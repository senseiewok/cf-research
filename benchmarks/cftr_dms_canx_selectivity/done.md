# Done file: Quantitative CFTR calnexin dependency (`CANX KO`) & single-corrector (`VX-661` vs. `VX-445`) selectivity benchmark (`Tedman et al. 2025`)

Pre-registered and committed before evaluating holdout outcomes or recording `ledger.json`.
Research benchmark only: plasma-membrane expression in HEK293T DMS assays is not channel conductance, clinical lung response, or medical advice.

## Success test

The benchmark is complete when all deterministic checks pass on the fixed protocol below without post-hoc retuning after holdout inspection:

1. **Data provenance, strict split reuse, & alignment policy:**
   - Source workbook: CC BY 4.0 `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` (`sha256`: `2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b`) in git-ignored `sources/downloads/tedman-2025/`, reading sheets `All_variants` (parental HEK293T) and `All_variants CANX KO` (calnexin-knockout HEK293T).
   - Strictly reuse the locked residue-clustered split from `benchmarks/cftr_dms_expression/split_manifest.json` (`split_manifest_sha256`: `e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543`): `167` `dev` variants across `123` unique codons (containing all `45` variants across the `30` exposed `PRED-0010` codons) and `67` `holdout` variants across `59` unique codons (`0` overlap with `PRED-0010`).
   - Require `234/234` identical `(row_num, variant, sub_domain)` alignment between `All_variants` and `All_variants CANX KO`; fail closed on any alignment mismatch.
   - Pre-registered positivity rule (no post-hoc imputation or floor tuning): require strictly positive finite values (`> 0`) in all 5 required Mut/WT relative intensity columns (`All_variants` Cols `F`, `J`, `N`, `R` and `All_variants CANX KO` Col `F`); fail closed if more than `5` variants must be excluded for non-positive or missing entries.

2. **Two independently gated scientific arms (Bonferroni family-wise control, no OR-gate):**
   - **Arm A (`canx_chaperone_coupling`):** Predict parental combination expression `y_combo = log2(Parental VX-661+VX-445)` using unpenalized 2-feature OLS `[1, log2(Parental_DMSO), log2(CANX_KO_DMSO)]` fitted on `dev` against Baseline `dmso_ols` (`[1, log2(Parental_DMSO)]`), evaluated by holdout MAE in `log2` units. (Reports linear-scale MAE and OLS coefficient SEs as context only.)
   - **Arm B (`single_corrector_selectivity`):** Define stabilized parental single-corrector selectivity strictly from raw Mut/WT relative intensities (ignoring author-derived summary columns `W`–`AC`):
     `S_eps = log2((VX445_rel + 0.05) / (VX661_rel + 0.05))`,
     with fixed `eps = 0.05` (Gemma 4 council stabilization against low-expression ratio blow-up; un-dampened `S_raw = log2(VX445_rel / VX661_rel)` and subdomain medians reported as non-gated context only). Fit partial-penalized Ridge on `dev` (`0` L2 penalty on intercept and `log2(Parental_DMSO)`; L2 penalty `alpha in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)` selected solely inside `dev` via 5-fold residue-clustered CV on `dev`-standardized CFTR sub-domain and variant-class indicators). Baseline for holdout comparison is **locked before seeing holdout** as whichever of `constant_dev_median(S_eps)` or `1D dmso_ols(S_eps)` achieves lower 5-fold residue-clustered `dev`-CV MAE.

3. **Independent decision thresholds per arm (`alpha_FWER <= 0.025` one-sided via 97.5% two-sided CI):**
   - Each arm is adjudicated independently (no cross-arm OR-gate).
   - For each arm evaluated once on the locked `67`-variant / `59`-codon holdout, compute `2,000`-resample residue-clustered bootstrap (`seed = 20261123`) of `delta_MAE = MAE_locked_baseline - MAE_candidate`, reporting both 95% (`[2.5%, 97.5%]`) and Bonferroni-adjusted 97.5% (`[1.25%, 98.75%]`) quantile intervals.
   - **Pre-registered pass criterion per arm (`supported`):** relative MAE reduction `(MAE_locked_baseline - MAE_candidate) / MAE_locked_baseline >= 0.10` **and** `bootstrap_ci975_low > 0.0`. Otherwise record that arm honestly as `refuted` (or `inconclusive` if `>= 10%` point reduction with `ci975_low <= 0.0`) in `results_2026-11-23.json`, `ledger.json`, and `failures.jsonl`.

4. **Independent Cholesky audit & repo checks:**
   - Separate code (`benchmarks/cftr_dms_canx_selectivity/independent_audit.py` importing nothing from `cftr_dms_canx_selectivity.py` or `tools/`) validates `done_sha256`, `split_manifest_sha256`, disjoint codons, `0` `PRED-0010` holdout overlap, and (when local `.xlsx` data is present) independently re-parses both workbook sheets and re-solves Arm A and Arm B via Cholesky $L L^T$ decomposition with float parity `tol <= 1e-6`.
   - Unit tests (`test_*.py`), `--check`, `independent_audit.py --check`, `check_claim_ledger.py --git`, `check_anti_hallucination.py`, `check_repo_organization.py`, `tools/sources` unit tests, `check_codeql_patterns.py`, and `check_claim_strength.py` all exit `0`.

## Positive control

Before reading holdout metrics:

1. Dual-sheet row/variant/sub-domain alignment across `All_variants` and `All_variants CANX KO` passes `234/234` with `<= 5` exclusions.
2. Development-set cross-background stability correlation passes `Pearson r(log2(Parental_DMSO), log2(CANX_KO_DMSO)) > 0.70`.
3. Development-set pharmacological directionality passes both monotonic chains: `median(VX-661+VX-445) > median(VX-445) > median(DMSO)` and `median(VX-661+VX-445) > median(VX-661) > median(DMSO)` in `All_variants`.
4. Shuffled-predictor negative controls inside `dev` (`seed = 20261124`: permuting `CANX_KO_DMSO` for Arm A and `sub_domain`/`variant_class` features for Arm B across `dev` rows) both achieve `< 10%` 5-fold `dev`-CV MAE reduction.
5. Synthetic unit tests check dual-sheet mismatch rejection, dampened `S_eps` stability, 59-cluster `97.5%` bootstrap bounds, shuffled negative-control rejection, and planted-tamper failure.

## Stopping rule

Single confirmatory evaluation on the locked `67`-variant / `59`-codon holdout. Zero retuning of pseudocount `eps = 0.05`, features, positivity rules, Ridge `alpha` grid, or split boundaries after observing holdout outputs. Any control failure stops execution before holdout reporting; any arm failing `< 10%` MAE reduction or `ci975_low <= 0.0` stays marked `refuted` in `results_2026-11-23.json`, `ledger.json`, and `failures.jsonl`.
