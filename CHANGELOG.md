# Changelog

All notable changes to bijli will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

- **fdtd**: CPML absorbing boundaries for `Fdtd3d` — `try_enable_cpml`, `try_enable_cpml_with`, `disable_cpml`, `cpml` (same validation as 2D: thickness ≥ 1 and `2·t + 3 ≤ min(nx, ny, nz)`; errors leave the boundaries unchanged). Roden–Gedney CPML on all six faces, including edges and corners: 12 Ψ fields (one per field component and normal derivative) stored only in the face slabs, with κ scaling and b/c profiles at the integer (E) / half-cell (H) Yee positions, applied as corrections after the main loops so the PEC path is unchanged. Measured against a 53³ PEC reference grid at 8 probes (faces, edge, corner) for a 10-cell layer: −80 dB (classic), −77 dB (CFS, κ_max = 3). Energy falls to 2e-8 of its peak once the pulse has left (PEC box: ≈ 0.1). Cost (new bench `fdtd/3d_step_30x30x30_cpml8`): 452 µs vs 298 µs for `fdtd/3d_step_30x30x30` (+51%, 15 s criterion run — the 8-cell slabs cover over half of a 30³ grid)

## [1.2.0] — 2026-11-09

Toolchain/dependency refresh plus a full physics and robustness audit of every module.
No public signatures were changed or removed; APIs whose correct answer did not fit
their signature are deprecated with a correctly-typed replacement. Several functions
now return corrected values or `Err` where they previously returned wrong numbers —
see **Breaking (behavioral)** for the migration notes.

### Breaking (behavioral — no signature changes)

Migration notes for consumers (kiran, joshua, prakash):

- **fdtd**: `ricker_wavelet` is now the true normalized Ricker wavelet (even, peak 1). The previous waveform (first derivative of a Gaussian, peak ≈ 0.858/τ) is available as `differentiated_gaussian` — switch to it to keep old source behavior.
- **fdtd**: `TfsfSource::incident_h` now returns H^inc at x = bx − ½ (the node the left-edge TFSF correction needs), was bx + ½.
- **fdtd**: `nf2ff_2d` returns non-negative far-field magnitudes √ρ|E_z| (was signed values) and returns `Err` for TE simulations, `box_margin = 0` or a non-positive frequency.
- **fdtd**: 2D/3D `set_material`/`set_permittivity`/`set_permeability` return `Err` when the region starts outside the grid or is inverted (was a silent no-op), and `Material::conductivity` is now honored (conductors were previously treated as vacuum).
- **fdtd**: `Fdtd2d::enable_cpml` now actually absorbs; invalid thickness (needs `≥ 1` and `2·t + 3 ≤ min(nx, ny)`) logs a warning and leaves the walls unchanged — use `try_enable_cpml` to get the error.
- **wave**: trig-free Fresnel functions (`*_direct`, `reflectance_unpolarized`, `transmittance_unpolarized`, `schlick_reflectance`, `reflectance_at_interface`) return `Err` for `cos θ_i` outside [0, 1] (beyond a 1e-6 clamp) and for non-positive/non-finite indices. Flip the normal (and swap n₁/n₂) for back-face hits.
- **wave**: `schlick_reflectance` values change for n₁ > n₂ (now uses cos θ_t; 1.0 under TIR).
- **wave**/**material**: material-based Fresnel functions use wave impedance and return `Err` for conductors (σ ≠ 0), ε_r ≤ 0 or μ_r ≤ 0.
- **scattering**: `mie` returns `Err` for Im(m) > 0 — bijli uses the e^{+iωt} convention, so absorbing media have Im(m) < 0; pass `m.conj()` for n + iκ data. Also rejects x > `MIE_MAX_SIZE_PARAMETER` (1e5) and |m|·x > `MIE_MAX_MX` (1e7).
- **beam**: `hermite_gaussian` returns the real amplitude envelope (identical at z = 0; differs for z ≠ 0). Use `hermite_gaussian_field` for the full complex field.
- **rf**: Touchstone 2-port files are read/written in the spec order `N11 N21 N12 N22` (S21/S12 were swapped); the parser rejects malformed data instead of silently shifting values; output numbers use scientific notation.
- **rf**: `single_stub_match` returns corrected values; a matched load now returns `(0.0, 0.25)` (λ/4 short-circuited stub ≡ no stub) instead of `(0.0, 0.0)`.
- **field**: `trace_field_line` returns the partial line (`Ok`) when `field_fn` reports `Singularity` or an infinite field (was `Err`).
- **field**: `image_charge_sphere` with the charge at the origin returns `InvalidParameter` (was `Singularity`).
- **field**: physical constants moved to CODATA 2018 (`EPSILON_0` −4.7e-10, `MU_0` +9.5e-11, `COULOMB_K` +2.9e-11 relative).

### Fixed

**fdtd**
- **fdtd**: `Fdtd2d` CPML was never applied — boundaries reflected 100% (`enable_cpml` built the layer but `step_tm`/`step_te` ignored it → Roden–Gedney CPML with graded σ/κ/CFS-α and strip-only Ψ fields at the correct integer/half-cell Yee positions, TM and TE; −88 dB reflection measured for a 10-cell layer)
- **fdtd**: `set_material` silently ignored `Material::conductivity` (no lossy update → per-cell C_a/C_b lossy E update in `Fdtd2d` and `Fdtd3d`, allocated only when σ > 0)
- **fdtd**: `Fdtd2d::set_material` left a partial update when μ_r was invalid (ε_r written before validation → all inputs and the region validated before any write, 2D and 3D)
- **fdtd**: stepping, `total_energy` and sources could panic out of bounds after the public `nx`/`ny`/`nz` fields or field vectors were mutated (→ layout checked each step; on mismatch, `tracing::error` and the step is skipped)
- **fdtd**: 2D wall comments said "Simple ABC" but walls reflect fully (→ documented as PEC in TM (E_z = 0), PMC in TE (H_z = 0); 3D documented as PEC)
- **fdtd**: `TfsfSource` auxiliary grid reflected the incident pulse at both ends (H copy acted as PEC, E copy as PMC → first-order Mur ABCs at both ends)
- **fdtd**: `TfsfSource` injected nothing when `boundary_x < 4` (source on node 0, overwritten by the boundary → clamped to an interior node)
- **fdtd**: `TfsfSource::new` panicked or diverged for `nx < 2`, `boundary_x ≥ nx`, `dt > dx/c` (→ sanitizes with warnings; `step`/`incident_*` never panic; `try_new` validates)
- **fdtd**: `Fdtd1d` Mur ABC coefficient was fixed at the vacuum value (≈ 11% energy reflection at ε_r = 4 → recomputed per end from the local ε_r·μ_r whenever materials change)
- **fdtd**: subpixel parallel/perpendicular formulas were swapped (E ∥ interface needs ⟨ε⟩, E ⊥ needs ⟨ε⁻¹⟩⁻¹ — Farjadpour et al. 2006); `subpixel_permittivity_parallel` now validates `fill_fraction ∈ [0, 1]`
- **fdtd**: `nf2ff_2d` dropped the imaginary part and the magnetic current M = −n̂×E, mispositioned H by half a cell and accepted TE simulations (→ Love-equivalence transform with J and M, complex accumulation, Yee-correct positions)
- **fdtd**: `fourth_order_derivative` wrote the right-edge value to index n−1 instead of n−2 (outputs now documented as derivatives at x_{i+½})
- **fdtd**: `find_eigenmodes` Q estimate scanned the whole spectrum (multi-mode Q off by ~100× → per-peak half-power walk with interpolation) and the frequency grid never reached `freq_max` (df = Δf/n → Δf/(n−1))
- **fdtd**: `ricker_wavelet` returned the first derivative of a Gaussian (→ true Ricker wavelet; see Breaking)
- **fdtd**: Drude ADE doc stated `β·E^{n+½}` but the coefficient is for the `β(E^{n+1}+E^n)` form (doc corrected; Drude/Lorentz/Debye now verified against analytic ε(ω) to < 1e-5)
- **fdtd**: `CwAccumulator::accumulate` panicked when the field was longer than the accumulator (→ common prefix with a one-time warning; `try_accumulate` errors)
- **fdtd**: `NonUniformGrid::courant_dt` returned ∞ for an empty axis and ≤ 0 for non-positive spacing (→ 0.0 with a warning; checked variants added)
- **fdtd**: NaN `dx`/ε/μ accepted by `Fdtd1d`/`Fdtd2d`/`Fdtd3d` constructors and setters (`<= 0.0` checks → finite and positive required; grid-size overflow returns `Err`); docs now say dt is *half* the Courant limit
- **fdtd**: τ = 0 produced NaN in the source functions (→ 0.0)

**compute**
- **compute**: `GpuFdtd2d::run` indexed out of bounds or wrapped to the next row for out-of-grid sources, and allocated a full host buffer every step (→ ignored like `Fdtd2d::add_source`; in-place point update)
- **compute**: `CpuBackend` panicked on buffer-length mismatches (reachable via `read_ez`/`read_hx`/`read_hy`) and underflowed for `ny = 0` (→ `InvalidParameter`); `GpuFdtd2d::new` rejects NaN `dx`
- **compute**: `ComputeBackend::fdtd_update_h` doc contradicted the reference implementation (axis/sign corrected); `WgpuBackend` documented as a non-functional scaffold

**field / charge / maxwell / relativity**
- **field**: `electric_quadrupole_potential` was 2× too small (divided by 4·4πε₀r⁵ → (1/4πε₀)·½ΣQᵢⱼxᵢxⱼ/r⁵ with Jackson's Qᵢⱼ = Σq(3xᵢxⱼ − r²δᵢⱼ), now documented)
- **field**: `electric_field_disk_axis` had the wrong sign/magnitude for z < 0 and returned exactly 0 for z ≫ R (z > 0-only formula, 1 − z/r cancellation → sgn(z)(σ/2ε₀)R²/(r(r+|z|)))
- **field**: `electric_field_ring_axis` raised false `Singularity` for atomic-scale rings (absolute 1e-30 floor on (z²+R²)^{3/2} → singular only for R = 0 at z = 0)
- **field**: `electric_field_inside_sphere` silently wrong for r outside [0, R]; `image_charge_sphere` accepted R ≤ 0/NaN; `image_charge_dielectric` accepted z = 0 and ε ≤ 0; Gauss-law fields accepted negative/NaN distances (→ `InvalidParameter`/`InvalidPermittivity`)
- **field**: `FieldVector::normalized` rejected legitimate SI magnitudes (< 1e-15, e.g. 1.6e-19 N) and returned a zero vector for components ≳ 1e154 (→ exact-zero check plus rescaled slow path)
- **field**: `trace_field_line` could overflow/OOM on huge `max_steps` and accepted a NaN step (→ capped initial capacity, NaN-safe validation)
- **field**: `multipole_moment` lost sin θ near the poles (acos→cos round trip → cos θ = z/r, sin θ = ρ/r); doc now states the unnormalized real moments it returns
- **field**: constants mixed CODATA sets, ε₀μ₀c² − 1 ≈ 3.8e-10 (→ CODATA 2018 ε₀, μ₀; `COULOMB_K` = 1/(4πε₀) computed at compile time; `bridge` constants updated to match)
- **relativity**: `solve_retarded_time` returned `NoConvergence` for t ≳ 1e-4 s and for some β ≥ 0.99 orbits (|f| < 1e-12 m stop was below c·ulp(t) rounding noise; unsafeguarded Newton → bracketed safeguarded Newton/bisection on τ = t − t_r with ulp-level stopping; t_r ≤ t guaranteed)
- **relativity**: `lienard_wiechert_fields` failed for fast orbits (hard-coded 50 Newton iterations → new solver, 200-step budget)
- **relativity**: `retarded_scalar_potential`/`retarded_vector_potential` docs claimed |r − r'(t_r)| though the 1/(1 − n̂·β) factor is absent (→ documented as the effective distance |R| − R·β; negative/NaN rejected; full LW functions added)
- **relativity**: `lorentz_transform_fields` ignored boosts with |v| < 1e-15 m/s (absolute shortcut → Jackson 11.149 form); NaN passed through `lorentz_factor`, `lorentz_transform_fields`, `relativistic_larmor_power`, `lienard_wiechert_e`
- **relativity**: `relativistic_larmor_power` doc said 6πε₀c (→ 6πε₀c³; code was correct)
- **charge**: `larmor_radius` rejected |q|B < 1e-30, e.g. electrons in B < 6e-12 T (product floor → separate validation); `cyclotron_frequency` accepted negative mass/|B|; `coulomb_potential_energy` accepted negative distance; `PointCharge::new` accepted NaN mass
- **maxwell**: `skin_depth` accepted two negative factors (sign of product only → each factor checked); all guards reject NaN

**wave / material / polarization / scattering / beam**
- **scattering**: `mie` was wrong for high-index or large spheres — x=20, m=4: Q_ext 2.0758 → 2.0890, g 0.4230 → 0.5009; x=50, m=3.5: 2.2122 → 2.0807; x=200, m=1.5: 2.0150 → 2.0921 (D_n(mx) downward recurrence started at n_max+15 regardless of |mx| → max(n_max, |mx|) + 4|mx|^{1/3} + 15; verified against a 40-digit reference)
- **scattering**: `mie` asymmetry g omitted the final (2n+1)/(n(n+1))·Re(aₙbₙ*) term
- **scattering**: `mie` Q_sca/Q_ext were 50× too large at x = 1e-8 (cancellation in ψ₁ = sin x/x − cos x → series for x < 0.1, Smith complex division)
- **scattering**: `mie` silently returned negative Q_abs for Im(m) > 0, panicked or allocated hugely for x = ∞ or x ≈ 1e10, and returned `Ok(NaN)` for NaN input (→ documented e^{+iωt} convention and validation; see Breaking)
- **scattering**: `rayleigh_angular_intensity` doc formula was 2× off and mixed polarizability units; a_n/b_n comment swapped (docs only)
- **beam**: `resonator_beam_waist` errored for the stable symmetric confocal cavity (0/0 at g₁ = g₂ = 0 → √(Lλ/2π))
- **beam**: `transverse_mode_spacing` was wrong for g₁, g₂ < 0 (dropped the ± of arccos(±√(g₁g₂)) → g = −0.5 gives 2·FSR/3, was FSR/3)
- **beam**: `hermite_gaussian` was neither |E| nor Re E away from the waist, e.g. HG₁₀ ≈ 0 everywhere at z = z_R (→ amplitude envelope; complex field via `hermite_gaussian_field`)
- **beam**: `AbcdMatrix::transform_beam` was wrong across media (assumed output medium = input → λ′ = λ·det(M), w₀′ = √(z_R′λ′/π))
- **wave**: `schlick_reflectance` underestimated R badly for n₁ > n₂ and stayed < 1 under total internal reflection (→ uses cos θ_t, returns 1.0 under TIR)
- **wave**: trig-free Fresnel family accepted any cos θ_i, e.g. `reflectance_unpolarized(1, 1.5, −1)` = 25 (→ range-checked; R, T guaranteed in [0, 1])
- **wave**: Fresnel `*_direct` returned `DivisionByZero` at grazing incidence on an index-matched interface (→ exactly 0 when n₁ = n₂)
- **wave**/**material**: `reflectance_at_interface`, `reflectance_normal_materials`, `material_reflection_coefficient` and `Material::reflectance_at` were wrong for magnetic media (ε = μ = 4 vs vacuum gave R = 0.36, true 0) and gave R = 0 for copper (n-only formula, σ ignored → impedance-form Fresnel; conductors rejected)
- **wave**: `rectangular_waveguide_cutoff(a, b, 0, 0, v)` returned `Ok(0)` and accepted any velocity (→ `Err`)
- **material**: `Material::impedance` accepted μ_r ≤ 0 (→ `InvalidPermeability`); `clausius_mossotti` returned negative ε_r as `Ok` when Nα/3ε₀ > 1 (→ `Err`)
- **polarization**: `JonesVector::normalized`/`overlap` and `StokesVector::degree_of_polarization` rejected valid tiny-amplitude states, and `normalized` returned the zero vector when the intensity overflowed (absolute 1e-30 thresholds → exact-zero checks with rescaling)

**rf / circuit / bridge / integration / ai**
- **rf**: `input_impedance_lossless` dropped the j and returned a wrong real value, e.g. 50 Ω instead of 40 − j30 Ω (→ new `input_impedance_lossless_complex`; the old function returns `Err` when Z_in is complex)
- **rf**: `single_stub_match` gave wrong positions/lengths (series-stub condition and a dimensionally wrong susceptance → Pozar §5.2 shunt short-stub solution, verified against Pozar Example 5.2)
- **rf**: Touchstone 2-port data was read and written row-major (v1.1 specifies `N11 N21 N12 N22`)
- **rf**: Touchstone files with ≥ 3 ports failed to parse (line-based parser → token stream across lines; writer wraps rows at 4 pairs per line)
- **rf**: Touchstone parser silently dropped bad tokens and short lines, honored later `#` lines, and ignored the parameter type and an invalid `R` (→ `Err`; inline `!` comments; 2-port noise block skipped; port count limited to 1..=256 with checked sizing so untrusted files cannot drive huge allocations)
- **rf**: `write_touchstone` lost precision for low frequencies and small S values and wrote `-inf` dB for |S| = 0 (`{:.6}` → shortest round-trip scientific notation; dB floor)
- **rf**: `SMatrix::get`/`set`/`cascade`, `s_to_z` and `write_touchstone` could panic on inconsistent public-field data (direct indexing → shape checks returning `Err`)
- **rf**: `cascade` and `s_to_z` silently used only `z0[0]` (unequal port reference impedances → `Err`)
- **rf**: `characteristic_impedance_lossy` returned NaN at DC; `gamma_from_z(−1)`/`z_from_gamma(1)` returned NaN at the Smith-chart poles (→ +∞ + 0j); `vswr`/`return_loss`/`insertion_loss` accepted NaN; `return_loss` accepted |Γ| > 1; `quarter_wave_transformer` accepted two negative impedances
- **circuit**: `quality_factor`, `damping_ratio` and `resonant_(angular_)frequency` returned `Ok(NaN)`/`Ok(0)` for unvalidated or negative L/C (→ R, L, C all validated; documented as series-RLC)
- **bridge**: `magnetic_flux_to_magnetostrictive_strain` kept growing beyond saturation (→ B/B_sat clamped); docs corrected for Joule-heating units, `em_energy_to_heating_rate` (takes W/m³, returns K/s), the eV↔nm constant (exact SI hc = 1239.841984 eV·nm), the BK7 validity range, the magnetostriction model and the piezo sign convention vs IEEE 176
- **integration/soorat**: `FieldSlice3D::from_fdtd_3d(sim, "Ex")` returned E_z data labelled "Ex" (case-sensitive match with fallback → case-insensitive; label names the data actually returned)
- **ai**: HTTP error statuses were treated as success, the configured `api_key` was never sent, `Debug` printed the key, and response bodies were unbounded (→ `error_for_status`, Bearer auth, redacted `Debug`, 1 MiB body cap)

**build**
- **features**: `charge`, `maxwell`, `material`, `relativity`, `wave`, `fdtd` and `compute` failed to compile on their own (`unresolved import crate::field` — missing inter-feature dependencies → declared in `Cargo.toml`; new `make features` target and CI job check every feature in isolation, including test code)

### Added

- **fdtd**: `CpmlParams` (thickness, grading order, σ_max/σ_opt ratio, κ_max, CFS α_max), `Cpml::with_params`, `Cpml::{sigma_max, kappa_max, alpha_max, grading_order}`
- **fdtd**: `Fdtd2d::{try_enable_cpml, try_enable_cpml_with, disable_cpml, cpml, set_conductivity}`, `Fdtd3d::{set_conductivity, set_permeability}`
- **fdtd**: `nf2ff_2d_phasor` — frequency-domain NF2FF from `CwAccumulator` phasors (validated against the exact Hankel-function line source to 0.4%)
- **fdtd**: `differentiated_gaussian`, `TfsfSource::try_new`, `try_new` for `LossyMaterial`/`DrudeMaterial`/`LorentzMaterial`/`DebyeMaterial`, `subpixel_permittivity_perpendicular_checked`, `NonUniformGrid::{try_uniform, validate, try_courant_dt}`, `CwAccumulator::{try_new, try_accumulate}`
- **compute**: `ComputeBackend::add_to_element` (provided method), `GpuFdtd2d::add_source`
- **field**: `spherical_multipole_moment` — normalized complex q_lm = Σ q rˡ Y*_lm (Jackson 4.3, Condon–Shortley phase), stable to l = 200
- **relativity**: `lienard_wiechert_scalar_potential`, `lienard_wiechert_vector_potential` (full LW potentials), `EmTensor::four_force` (Lorentz four-force)
- **beam**: `hermite_gaussian_field` (complex HG field, e^{−iωt}), `AbcdMatrix::transform_beam_with_waist_position`
- **scattering**: `MIE_MAX_SIZE_PARAMETER`, `MIE_MAX_MX`
- **material**: `Material::try_reflectance_at`
- **rf**: `input_impedance_lossless_complex`, `single_stub_match_complex` (both solutions, complex loads), `try_write_touchstone`
- **integration/soorat**: `FieldSlice2D::try_from_fdtd_2d`, `FieldSlice3D::try_from_fdtd_3d`, `VectorFieldSample::try_from_components`
- **benches**: 15 new benchmarks — CPML/lossy 2D step, 3D step, NF2FF, eigenmodes, Mie x=50, HG complex field, material Fresnel, complex input impedance, stub match, Touchstone parse/write, retarded-time solve at t = 1 s, LW fields on a β = 0.99 orbit, LW scalar potential

### Changed

- **toolchain**: MSRV 1.89 → 1.99 (latest stable; edition 2024 remains the newest edition)
- **relativity**: `EmTensor`/module docs state that stored F^μν uses the (+,−,−,−) Jackson/Griffiths layout, i.e. −F^μν under `FourVector`'s (−,+,+,+) metric (component values unchanged)
- **relativity**: `solve_retarded_time` no longer uses `hisab::num::newton_raphson` (supersedes the 1.1.0 entry); `max_iterations` counts refinement steps
- **polarization**: `ComplexExt::norm_sq`/`arg` and `complex_from_polar` delegate to hisab 1.4's inherent `Complex` methods (identical results)
- **polarization**/**scattering**/**beam**: time conventions documented (polarization/beam e^{−iωt}, scattering e^{+iωt}), plus waveplate fast/slow phases and `MuellerMatrix::rotation(θ)` = Mueller form of `JonesMatrix::rotator(−θ)`
- **rf**: module docs list the conversions that actually exist (2-port cascade via internal T-parameters, 2-port S→Z; the 1.0.0 entry's "N-port S↔Z↔T" overstated this)
- **bridge**: `em_energy_to_heating_rate` parameter renamed `power_density_w_m3`
- **ci**: actions/checkout v4 → v7, codecov-action v4 → v7, action-gh-release v2 → v3, MSRV job → 1.99; new Feature Isolation job
- **deny**: removed the stale `RUSTSEC-2024-0436` (paste) ignore and the deprecated `GPL-3.0` license allowance — neither is reachable with the updated dependency tree

### Deprecated

- **rf**: `input_impedance_lossless` → use `input_impedance_lossless_complex`
- **rf**: `write_touchstone` → use `try_write_touchstone`
- **material**: `Material::reflectance_at` → use `try_reflectance_at`

### Dependencies

- hisab 1.3 → 1.4
- wgpu 25 → 30 (optional, `compute-wgpu`; drops the unmaintained `paste` transitive dependency)
- Minimum versions raised to current: serde 1.0.229, thiserror 2.0.21, tracing 0.1.44, reqwest 0.13.5, tokio 1.54, serde_json 1.0.154, tracing-subscriber 0.3.23, criterion 0.8.2

### Performance

Measured as a quiet back-to-back A/B (pre-audit tree with updated deps vs. this release, same machine, criterion medians; only changes with non-overlapping confidence intervals are listed as wins/regressions):

- `beam/hermite_gaussian_hg33`: 22.70 ns → 5.17 ns (−77%) — envelope no longer evaluates atan/cos
- `fdtd/2d_step_100x100_tm`: 15.59 µs → 12.94 µs (−17%) — row-slice update loops
- `fdtd/2d_run_50_steps_80x80`: 597.39 µs → 468.99 µs (−21%)
- `scattering/mie_x10_real`: 996.02 ns → 892.22 ns (−10%) — single-division reciprocal in the D_n recurrence, 1/m hoisted (despite the more accurate recurrence start)
- `scattering/mie_x1_real`: 484.15 ns → 465.16 ns (parity, within noise)
- `fdtd/step_1000_cells`, `step_100_cells`, `run_100_steps_500_cells`, `wave/reflectance_unpolarized`: parity (within noise) despite added validation
- **Regression** — `rf/s_cascade`: 77.24 ns → 94.10 ns (+22%) — cost of the new shape and reference-impedance consistency checks that replace panics/silent wrong results; reciprocals switched to a single-division fast path to limit it
- New paths: `fdtd/2d_step_100x100_tm_cpml10` ≈ 1.8× the PEC step (CPML strips cover 36% of the 100×100 grid); `fdtd/2d_step_100x100_tm_lossy` ≈ 1.2× (per-cell C_a, allocated only when σ > 0)

## [1.1.0] — 2026-03-29

### Added
- **bridge** — cross-crate primitive-value bridges for ushma (Joule heating, EM energy to heating rate), prakash (Sellmeier refractive index, E-field to optical intensity), tanmatra (energy level to wavelength, nuclear charge potential), dravya (piezoelectric stress, magnetostrictive strain)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `FieldSlice2D`/`FieldSlice3D` (FDTD heatmaps with `from_fdtd_2d`/`from_fdtd_3d`), `FieldLineVisualization` (polylines with magnitudes), `ChargeVisualization` (point charges with `from_point_charges`), `RadiationPattern` (far-field with `from_nf2ff`), `VectorFieldSample` (arrow/streamline grids)

### Changed

- **deps**: iri-string 0.7.11 → 0.7.12
- **license**: GPL-3.0 → GPL-3.0-only (no "or later" clause)
- **deps**: hisab 1.1.0 → 1.3.0
- **polarization**: `Complex` type now re-exported from `hisab::num::Complex`; added `ComplexExt` trait and `complex_real`/`complex_zero`/`complex_from_polar` free functions
- **relativity**: `solve_retarded_time` now uses `hisab::num::newton_raphson` for root-finding

### Removed

- **deps**: Removed unused `impetus` dependency
- **polarization**: Removed custom `Complex` implementation (~185 lines) in favor of hisab's

## [1.0.0] — 2026-03-25

### Added

**New modules:**
- **polarization**: Jones vectors/matrices (coherent light), Stokes parameters, Mueller matrices (partially polarized light), optical elements (polarizer, waveplate, retarder, rotator, depolarizer), Complex number type
- **scattering**: Mie scattering (Bohren-Huffman algorithm, a_n/b_n coefficients, Q_sca/Q_ext/Q_abs, asymmetry parameter g), Rayleigh scattering (λ⁻⁴ cross-section, phase function)
- **beam**: Gaussian beam parameter q(z), ABCD ray transfer matrices, resonator stability (g-parameters), Hermite-Gaussian and Laguerre-Gaussian mode profiles
- **rf**: Transmission lines (lossless/lossy Z₀, propagation constant), VSWR, return/insertion loss, Smith chart, matching networks (quarter-wave, single-stub), S-parameters (N-port, cascade, S↔Z↔T conversions), Touchstone I/O, antenna arrays (linear, planar, circular)
- **compute**: `ComputeBackend` trait (vendor-agnostic GPU/CPU), `CpuBackend` reference impl, `GpuFdtd2d<B>` generic solver, `SooratBackend`/`WgpuBackend` scaffolds

**Existing module additions:**
- **wave**: Trig-free Fresnel equations (`fresnel_rs_direct`, `fresnel_rp_direct`, `snell_cos_theta_t`), Schlick approximation, material-based Fresnel interface
- **fdtd**: 2D FDTD (TM/TE modes), 3D FDTD (full vector Yee grid, 6 field components), CPML absorbing boundaries, source functions (Gaussian pulse, Ricker wavelet, modulated Gaussian, TFSF), dispersive materials (Drude/Lorentz/Debye with ADE coefficients, gold/silver/water presets), NF2FF transform, subpixel smoothing, `NonUniformGrid`, 4th-order spatial stencil, `CwAccumulator`, `find_eigenmodes`
- **field**: Multipole expansion (quadrupole potential, spherical harmonic moments, associated Legendre polynomials), method of images (conducting plane, grounded sphere, dielectric half-space), Green's functions (scalar complex + static), Maxwell stress tensor, Poynting flux surface integration
- **relativity**: Retarded time solver (Newton-Raphson), convenience `lienard_wiechert_fields()` combining retarded time + E + B
- **charge**: `PointCharge` convenience methods — `electric_field_at`, `electric_potential_at`, `magnetic_field_at`, `lorentz_force_in`
- **material**: Unified `Material` struct — vacuum/dielectric/conductor constructors, refractive index, impedance, wave speed, reflectance

### Changed

- **circuit**: `resistance_series`, `capacitance_parallel`, `inductance_series` now return `Result` and reject empty slices
- **fdtd**: Mur ABC coefficient precomputed in constructor; added tracing instrumentation
- **fdtd**: `subpixel_permittivity_parallel` now returns `Result` (validates permittivity > 0)
- **fdtd**: `Fdtd2d::set_material` and `Fdtd3d::set_material` accept unified `Material` struct
- **rf**: `input_impedance_lossless` now returns `Result` (catches zero denominator)
- **rf**: `normalize_impedance` now returns `Result` (validates Z₀ > 0)
- **field**: Added `Default`, `From<[f64;3]>`, `Into<[f64;3]>` for `FieldVector`; exposed `SPEED_OF_LIGHT_SQ`
- **field**: Removed redundant `FieldVector::add()` method (use `+` operator)
- **polarization**: Added `Default` for `Complex`, `JonesVector`, `StokesVector`; `From<f64>` and `Div` for `Complex`
- **polarization**: Removed redundant `StokesVector::add()` method (use `+` operator)

### Dependencies

- hisab 0.22.4 → 1.1.0
- impetus 0.23.3 → 1.0.0
- reqwest 0.12 → 0.13
- criterion 0.5 → 0.8
- Added: wgpu 25 (optional, `compute-wgpu` feature)
- Removed: soorat dependency

## [0.24.3] — 2026-03-24

### Added

- **field**: Electric/magnetic dipole fields (full 3D), dipole potential
- **field**: Gauss's law applications — infinite plane, charged sphere (inside/outside), infinite line, cylinder (inside/outside)
- **field**: Charge distributions — ring (axial), disk (axial)
- **field**: Field line tracing with Euler integration
- **field**: `std::ops` traits (`Add`, `Sub`, `Mul<f64>`, `Neg`) and `PartialEq` on `FieldVector`
- **relativity**: New module — Lorentz factor, beta, four-vectors, electromagnetic field tensor F^μν, Lorentz transformation of fields, retarded potentials, Liénard-Wiechert fields, Larmor radiation power
- **wave**: Snell's law, Fresnel equations, waveguides (rectangular/cylindrical), antenna patterns, Friis transmission
- **fdtd**: New module — 1D FDTD Yee algorithm with Mur ABC
- **circuit**: New module — Ohm's law, series/parallel combinations, RC/RL/RLC circuits
- **material**: New module — dielectric/magnetic materials, Clausius-Mossotti, Curie's law

### Changed

- **charge**: `PointCharge::new` now validates mass > 0
- **maxwell**: `refractive_index`, `wavelength`, `frequency` now validate inputs

### Performance

- `inv_r` multiply pattern — 10-29% faster field/charge computations
- `superposition_10`: 64ns → 41ns (-36%)

## [0.1.0] — 2026-03-24

### Added

- **field**: Electric/magnetic vector fields, point charge fields, potentials, superposition, energy density, constants (ε₀, μ₀, c, k)
- **maxwell**: Gauss's laws, displacement current, wave speed, refractive index, impedance, wavelength/frequency, skin depth
- **charge**: Point charges (electron, proton), Coulomb force, Lorentz force, dipole moment, potential energy, cyclotron frequency, Larmor radius
- **wave**: Poynting vector, plane wave intensity, radiation pressure, E/B amplitude conversion, wave number, momentum density
- **error**: BijliError with domain-specific variants
- **ai**: Daimon/hoosh client integration (feature-gated)
- **logging**: Structured logging via BIJLI_LOG (feature-gated)
