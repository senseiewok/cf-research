# `benchmarks/mucus_transport/` — Correctness, Uncertainty Split & Identifiability Audit (`Items #6 & #8`)

Independent narrow-correctness, dimensional, source-split uncertainty, and algebraic parameter identifiability audit accompanying [tools/discovery/mucus_barrier_screener.py](../../tools/discovery/mucus_barrier_screener.py).

## Scope of This Audit

- **What this audit checks:**
  1. **Arithmetic parity & SI dimensions (`audit_dimensional_consistency`, `audit_parity_and_invariants`):** Re-derives the Renkin-Flory-Stokes-Einstein screener equations independently, checks $(M, L, T, \Theta)$ exponents (`D_w, D_eff -> L^2 T^-1`; `transit_time -> T`) and the $\mu\text{m} \to \text{m}$ scaling factor, pins the documented cationic-stricter-than-anionic charge asymmetry (`+zeta` vs `-zeta`), checks smooth-domain monotonicity and active input dependence, and flags the `C0` step jumps at $R_F / D \in \{1.0, 1.5\}$ plus the `C1` slope kink at $|\zeta_{\text{eff}}| = 5\text{ mV}$.
  2. **Five planted defects (`run_planted_defect_checks`):** Exercises detection of missing $\mu\text{m} \to \text{m}$ unit conversion ($10^{12}$ magnitude error), electrostatic exponent sign flip, missing steric exclusion clamp, disconnected pore diameter input, and swapped charge-sign asymmetry.
  3. **Source-split uncertainty (`run_source_split_uncertainty_audit`):** Keeps literature bounds (`pore_diameter_nm`, `mucus_thickness_um`, `mcc_window_s`) strictly separated from uncalibrated heuristic parameters (`interstitial_viscosity_ratio`, `charge_shield_*`, `electrostatic_alpha_*`) across `measured_only`, `heuristic_only`, and `joint_measured_plus_heuristic` grids without merging them into one unlabeled interval. Handles steric boundary singularities (`transit_time_is_finite: false`) cleanly without `Inf`/`NaN` serialization.
  4. **Synthetic known-answer parameter identifiability (`identifiability_proofs`):** Tests structurally that single-size neutral particle diffusivity confounds `(d_pore, eta_ratio)`, that two distinct neutral probes $(0 < d_1 < d_2 < d_{\text{pore}})$ algebraically recover both `(d_pore, eta_ratio)`, and that single charged probes confound `(charge_shield_factor, electrostatic_alpha)`.
- **What this audit does not claim:**
  - Arithmetic consistency and `SYNTHETIC_KNOWN_ANSWER` algebraic identifiability are **not** empirical validation in human CF sputum. No public cell-free or patient CF sputum particle-tracking dataset has been admitted yet (see [proposals/2026-11-22-mucus-mpt-data-feasibility.md](../../proposals/2026-11-22-mucus-mpt-data-feasibility.md)).
  - $f_{\text{steric}} = (1 - d_p / d_{\text{pore}})^2$ captures geometric steric exclusion without full Renkin-Hille hydrodynamic pore drag terms.

## Run Checks

```bash
python benchmarks/mucus_transport/audit_mucus_barrier.py --check
python benchmarks/mucus_transport/audit_mucus_barrier.py --check --json
python -m unittest discover -s benchmarks/mucus_transport -p 'test_*.py' -v
```
