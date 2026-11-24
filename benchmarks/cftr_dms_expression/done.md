# Done file: Quantitative CFTR plasma-membrane expression benchmark (`Council Item #1` / Tedman et al. 2025)

Pre-registered and committed alone before extracting holdout outcomes or recording `ledger.json`.
Research benchmark only: plasma-membrane expression in HEK293T DMS assays is not channel conductance, clinical lung response, or medical advice.

## Success test

The benchmark is complete when all deterministic checks pass on the fixed protocol below without post-hoc retuning:

1. **Data provenance & scope:** Uses maintainer-ratified CC BY 4.0 files from Tedman et al. (*eLife* 2025, `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`) placed locally in git-ignored `sources/downloads/tedman-2025/`:
   - `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` (`sha256`: `2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b`), sheet `All_variants` (parental cells).
   - `AGM ALT Library Comparison.xlsx` (`sha256`: `6673ab2a480200e7d36133b98dfbeb4864d1a41aa22cf4fdddea59eda8dbe178`), sheet `Sheet1`.
2. **Multiplicative log2 scale & zero imputation:**
   - Primary target: `y_combo = log2(VX-661+VX-445 Average Relative Intensity (Mut/WT))`.
   - Baseline untreated feature: `x_dmso = log2(DMSO Average Relative Intensity (Mut/WT))`.
   - Any row with missing, non-numeric, `<= 0` ratio for `DMSO` or `VX-661+VX-445`, or unparseable codon position is excluded cleanly with explicit counts and IDs logged (zero imputation).
   - Primary evaluation metric: mean absolute error (MAE) in `log2` units (fold-change bits), with linear-scale MAE reported alongside as secondary context only (never used for the decision gate).
3. **Zero-leakage residue-clustered split:**
   - Group all variants strictly by CFTR `codon_position` parsed from `Protein distinction` across missense (`R3W`), nonsense (`E60X`, `W1282*`), in-frame indel (`F508del`, `I507del`), and other classes.
   - Force every codon overlapping the 30 previously exposed `PRED-0010` benchmark variants (`benchmarks/theratype/theratype_open_annotations_2026-11-21.json`) into `dev` (`0` overlap with `holdout`).
   - Assign remaining unexposed residue clusters deterministically by codon hash: `int(sha256(f"tedman2025-item1:{codon}").encode("utf-8")[:8], 16) % 100 < 60` maps to `dev`, else `holdout`.
   - Require `holdout` to contain at least `40` unique unexposed residue clusters before proceeding.
4. **Fixed baselines vs. partial-penalized Ridge candidate (fitted on `dev` only):**
   - **Baseline 1 (`constant_dev_median`):** Median of `y_combo` on `dev`.
   - **Baseline 2 (`dmso_ols`):** Unpenalized ordinary least squares `y_combo = b0 + b1 * x_dmso` fitted on `dev`.
   - **Candidate (`dmso_subdomain_class_ridge`):** Partial-penalized Ridge regression `y_combo = b0 + b1 * x_dmso + sum_j w_j * z_j`, where `z_j` are `dev`-standardized indicators for CFTR sub-domain (`Lasso`, `TMD1`, `NBD1`, `R Domain`, `TMD2`, `NBD2`) and variant class (`missense`, `nonsense`, `inframe_indel`, `other`). Intercept `b0` and `x_dmso` slope `b1` carry `0` L2 penalty; L2 penalty `alpha` applies only to `w_j` and is chosen solely within `dev` via 5-fold residue-clustered CV over `(0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)`. Log the full CV path and compare candidate `b1` against Baseline 2 `b1`.
5. **One-shot holdout decision gate & independent verification:**
   - Evaluate Baseline 1, Baseline 2, and Candidate **once** on the locked `holdout`. Let `stronger_baseline` be the lower holdout MAE between Baseline 1 and Baseline 2 (also recording which won in dev CV).
   - Compute 2,000-resample seeded (`seed=20261123`) residue-clustered bootstrap 95% confidence interval on `delta_MAE = MAE_stronger_baseline - MAE_candidate`, explicitly reporting `n_holdout_clusters`.
   - **Pre-registered pass criterion (`supported`):** relative MAE reduction `(MAE_stronger_baseline - MAE_candidate) / MAE_stronger_baseline >= 0.10` **and** `bootstrap_ci95_low > 0.0`. Otherwise record `"refuted"` (or `"inconclusive"` if `CI` straddles zero despite point reduction) honestly in `ledger.json` and `failures.jsonl`.
   - Separate code (`benchmarks/cftr_dms_expression/independent_audit.py` importing nothing from `cftr_dms_benchmark.py` or `tools/`) must confirm `0` count discrepancy and `< 1e-6` float parity on locally parsed `.xlsx` data (and split/result invariants in CI when git-ignored downloads are absent).
   - `python -m unittest discover -s benchmarks/cftr_dms_expression -p "test_*.py" -v`, `check_anti_hallucination.py`, `check_repo_organization.py`, `check_codeql_patterns.py`, `check_claim_strength.py`, and `check_claim_ledger.py --git benchmarks/cftr_dms_expression/ledger.json` all exit `0`.

## Positive control

Before touching holdout evaluation:

1. Synthetic sparse `.xlsx` contract test recovers shared strings (`t="s"`), inline strings (`t="inlineStr"`), scientific-notation numbers, blank cells, and absent sparse `<c>` elements (`r="D3"`).
2. Partial-penalized Ridge synthetic control recovers known `b1` slope (`tol < 1e-6`) and planted feature weights without shrinking the `x_dmso` slope.
3. Cross-study control on `AGM ALT Library Comparison.xlsx` (`N = 103` rows comparing McKee et al. 2023 [AGM] vs. Tedman et al. 2025 [ALT] DMS intensities) reproduces strong positive Pearson `r > 0.85` and Spearman `rho > 0.85`.
4. Development-set directionality check (`dev` rows only): median untreated `DMSO` Mut/WT `< 1.0`, and median `VX-661+VX-445` Mut/WT `> DMSO` median.

## Stopping rule

Single confirmatory pass on the locked `holdout` residue clusters. Zero retuning of features, log transforms, codon exclusions, regularization grid, or split boundaries after inspecting holdout outcomes. If `AGM ALT` positive control fails (`r <= 0.85`), fewer than `40` unexposed holdout residue clusters remain, or the gate fails (`< 10%` MAE improvement or bootstrap `95% CI` lower bound `<= 0`), stop and publish the negative or blocked record cleanly in `results_2026-11-23.json`, `ledger.json`, and `failures.jsonl`.
