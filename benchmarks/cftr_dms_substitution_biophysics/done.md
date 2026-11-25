# Done file: Substitution Biophysics + UniProt `P13569` Topology DMS Rescue Benchmark (`cftr_dms_substitution_biophysics`)

Pre-registered and committed as a standalone commit before running outcome evaluation, reading holdout intensities, or recording `ledger.json`.
Research benchmark only: plasma-membrane expression in HEK293T DMS assays is not channel conductance, clinical lung response, or medical advice.

## Success test

The benchmark is complete when all deterministic checks pass on the fixed protocol below without post-hoc retuning:

1. **Data provenance, locked split, and offline UniProt `P13569` annotation:**
   - **DMS workbook:** Maintainer-ratified CC BY 4.0 file `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` (`sha256`: `2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b`, Tedman et al., *eLife* 2025, `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`) located in git-ignored `sources/downloads/tedman-2025/` (resolved across main checkout or `.worktrees/cf-research/...`, with graceful offline CI mode when absent). Parse sheet `All_variants` via sparse XML coordinates: `y_combo = log2(VX-661+VX-445 Mut/WT)` (Col `R`) and `x_dmso = log2(DMSO Mut/WT)` (Col `F`), requiring strictly positive finite values (`> 0`) and zero outcome imputation (`234/234` rows preserved).
   - **Locked residue-clustered split:** Reuse `benchmarks/cftr_dms_expression/split_manifest.json` (`sha256`: `e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543`): `167` `dev` variants across `123` unique codons (including all `45` variants across the `30` exposed `PRED-0010` codons) and `67` `holdout` variants across `59` unique codons (`0` overlap with `dev` or `PRED-0010`).
   - **Tracked offline UniProt `P13569` JSON (`tools/cftr-model/data/p13569.json`, CC BY 4.0):** `sha256`: `c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893` (verified against `tools/cftr-model/data/ATTRIBUTION.md`), `1,480` amino acids, parsing `Transmembrane`, `Helix`, `Beta strand`, and `Binding site` feature intervals.

2. **Variant parsing, experimental-WT rule, and 7-feature substitution biophysics vector:**
   - Handle all `234` variants:
     - **Single missense (`R3W`):** one substitution component `(wt_exp, pos, mut)`.
     - **Compound missense (`Row 77 Q359K/T360K`):** parse all `/`-separated components; sum continuous substitution deltas across components and take the logical-OR of binary topology indicators across component positions.
     - **Synonymous (`E528E`, `T854T`, `Q966=`, `L1156=`):** set all `7` features to `0.0`.
     - **In-frame indels (`L138ins`, `F312DEL`, `I507del`, `F508del`):** `is_inframe_indel = 1.0`, `is_nonsense = 0.0`, `helix_breaker_in_helix = 0.0`, `tm_polarity_shock = 0.0`, and continuous substitution scales set to missing (`None`) for in-fold `dev`-mean imputation.
     - **Nonsense (`E60X`..`Q1313X`):** `is_nonsense = 1.0`, `is_inframe_indel = 0.0`, `helix_breaker_in_helix = 0.0`, `tm_polarity_shock = 0.0`, and continuous substitution scales set to missing (`None`) for in-fold `dev`-mean imputation.
   - **Experimental-WT rule:** Compute substitution deltas strictly from the variant's experimental WT amino acid (`M -> V` for `Row 89 M470V`), rather than canonical `P13569`'s `V470`.
   - **Fixed 7-feature vector:**
     1) `delta_blosum62 = BLOSUM62(wt_exp, mut) - BLOSUM62(wt_exp, wt_exp)` using the published Henikoff & Henikoff (1992, *PNAS* 89:10915) `20x20` BLOSUM62 matrix (`W,W=11`, `C,C=9`, etc.)
     2) `delta_kd_hydropathy = KD(mut) - KD(wt_exp)` using canonical signed Kyte & Doolittle (1982, *J. Mol. Biol.* 157:105) hydropathy values (`I=+4.5` .. `R=-4.5`)
     3) `abs_delta_charge = |charge(mut) - charge(wt_exp)|` (`K,R=+1`, `D,E=-1`, others `0`)
     4) `helix_breaker_in_helix = 1.0` if `P13569` marks codon in `Helix` or `Transmembrane` AND missense/compound-missense has component with `wt_exp in {'P','G'}` or `mut in {'P','G'}`, else `0.0`
     5) `tm_polarity_shock = 1.0` if `P13569` marks codon in `Transmembrane` AND missense/compound-missense has component with `|delta_kd| >= 3.0` or `abs_delta_charge >= 1.0`, else `0.0`
     6) `is_nonsense` (`1.0` if nonsense, else `0.0`)
     7) `is_inframe_indel` (`1.0` if in-frame indel, else `0.0`)

3. **Baselines, partial-penalized Ridge candidate, and zero-leakage imputation:**
   - **Baselines (fitted on `dev`):** `constant_dev_median` and `dmso_ols` (`y_combo = b0 + b1 * x_dmso`), locking `stronger_baseline` by lower 5-fold codon-clustered `dev`-CV MAE.
   - **Candidate (`dmso_subst_biophysics_ridge`):** Partial-penalized Ridge regression (`0` L2 penalty on intercept `b0` and `x_dmso` slope `b1`; L2 penalty `alpha in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)` chosen solely via 5-fold `dev` codon-clustered CV on `dev`-standardized features). Missing continuous substitution scales (`None` on nonsense and in-frame indel rows) are imputed strictly to the training-fold mean inside each `dev` CV fold, and to the full `dev` mean for `holdout` (zero cross-split or cross-fold leakage).

4. **One-shot holdout decision gate, context diagnostics, & independent verification:**
   - Evaluate `constant_dev_median`, `dmso_ols`, and `dmso_subst_biophysics_ridge` once on the locked `67`-variant (`59`-codon) `holdout`.
   - Compute a `10,000`-resample seeded (`seed=20261124`) `59`-codon cluster-bootstrap 95% confidence interval on `delta_MAE = MAE_stronger_baseline - MAE_candidate`.
   - **Pre-registered pass criterion (`supported`):** `relative_mae_reduction >= 0.10` AND `bootstrap_ci95_low > 0.0`. Otherwise record `"refuted"` (or `"inconclusive"` if `relative_mae_reduction >= 0.10` while `bootstrap_ci95_low <= 0.0`) honestly in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`.
   - Report non-gated context explicitly: fitted `dev` weights (including `is_inframe_indel`), all `dev` feature VIFs (flagging if any `VIF > 10`), per-`sub_domain` holdout breakdown, and comparison vs `#396` (`+3.70%`) and `#401` (`+0.1261%`).
   - Independent audit (`benchmarks/cftr_dms_substitution_biophysics/independent_audit.py`, zero imports from `cftr_dms_substitution_biophysics.py` or `tools/`) independently parses `p13569.json` and `.xlsx` (when present locally; CI checks JSON/BLOSUM62/P13569/split invariants when absent), verifies `done_sha256`, `split_manifest_sha256`, and `p13569_sha256`, and solves partial Ridge via Cholesky `L L^T` (`<= 1e-6` float parity plus propagated rounding bound on printed relative reduction).
   - `unittest -v`, `--check`, `independent_audit.py --check`, `check_claim_ledger.py --git`, `check_claim_strength.py`, `check_codeql_patterns.py`, `check_anti_hallucination.py`, `check_repo_organization.py`, and `tools/sources` tests all exit `0`.

## Positive control

Before evaluating holdout outcomes, all three controls must pass (fail closed before reading holdout outcomes if any fails):

1. **`P13569` WT sequence parity across all 234 variants (including `Q359K/T360K`):** Across all components of all `234` variants, parsed experimental WT residues match `P13569` `sequence.value[codon - 1]` with zero unexpected mismatches anywhere (`0` in `holdout`) and only two documented `dev` exceptions: `(89, 470, 'M470V')` (`M` vs `V470`) and `(217, 966, 'Q966=')` (`Q` vs `T966`).
2. **`#396` baseline parity (`<= 1e-6`):** `dmso_ols` `dev` fit (`b0`, `b1`, `constant_dev_median`, 5-fold `dev`-CV MAE) and holdout MAE (`0.451487`) match `benchmarks/cftr_dms_expression/results_2026-11-23.json` within `1e-6`.
3. **Null biophysics control (`dev`-only, `seed=20261124`):** Permuting the 7-feature rows across `dev` variants yields `< 10%` 5-fold `dev`-CV MAE reduction relative to `stronger_baseline`.

## Stopping rule

Single-shot evaluation on the locked `67`-variant (`59`-codon) `holdout` after all pre-holdout controls pass. Zero retuning of variant parsing rules, BLOSUM62/Kyte-Doolittle/charge tables, topology thresholds, imputation rules, Ridge `alpha` grid, or split boundaries after reading holdout outputs. Any control failure stops execution before reading holdout; any miss against `>= 10%` MAE reduction or `bootstrap_ci95_low > 0.0` is recorded as `refuted` (or `inconclusive`) in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`.
