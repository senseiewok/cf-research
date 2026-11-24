# CFTR Calnexin (`CANX KO`) Chaperone Coupling & Single-Corrector (`VX-661` vs. `VX-445`) Selectivity Benchmark

This directory implements a pre-registered follow-up empirical benchmark on the **Tedman et al. 2025** CFTR deep mutational scanning dataset (*eLife*, DOI `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`, CC BY 4.0).

**Epistemic & clinical scope boundary:** Measurements reflect **plasma-membrane immunostaining expression (`Mut/WT` relative intensity)** in parental and calnexin-knockout (`CANX KO`) HEK293T cells. Plasma-membrane abundance is **not** CFTR chloride channel conductance (`PRED-0010`), patient airway physiology, or clinical treatment guidance.

## Pre-Registered Protocol & Multi-Family Council Guards

`done.md` was committed prior to holdout evaluation (`done_sha256: 9d701976cea245957e18b3feb532d6ad7b01338928f40397bb5472116aa81187`), reusing the locked residue-clustered split from `benchmarks/cftr_dms_expression/split_manifest.json` (`167` dev variants / `123` codons including all `45` variants across the `30` `PRED-0010` codons; `67` holdout variants / `59` unique codons; `0` overlap):

- **Workbook finding documented in `done.md`:** In the published Tedman et al. 2025 workbook, Parental `All_variants` `Col J` (`VX-661 Mut/WT Relative Intensity`) is blank (`0/234`), whereas `All_variants CANX KO` populates all four drug arms (`Cols F, J, N, R`, `234/234`).
- **Pre-build council protections (`ornith-1.5:35b` + `gemma4:31b-it-q4_K_M`):**
  1. Two independent gates with Bonferroni-adjusted two-sided **97.5%** residue-clustered bootstrap CIs (`alpha_FWER <= 0.025` one-sided, `2,000` resamples, `seed = 20261123`), no OR-gate across arms.
  2. Dampened selectivity target `S_eps_canx_ko = log2((CANX_KO_VX445_rel + 0.05) / (CANX_KO_VX661_rel + 0.05))` (`eps = 0.05`) to prevent low-expression ratio spikes.
  3. Dev-CV locking of the Arm B baseline before holdout inspection (`constant_dev_median` won dev CV over `canx_ko_dmso_ols`).
  4. Shuffled-predictor negative controls inside 5-fold dev CV (`seed = 20261124`).

## Observed Results (Single-Shot Holdout Evaluation)

- **Pre-holdout controls (all passed):** `234/234` dual-sheet alignment (`0` excluded), `dev` `log2(Parental DMSO)` vs. `log2(CANX KO DMSO)` Pearson `r = 0.992301`, monotonic `dev` medians in both Parental (`0.408232 < 1.462794 < 3.120957`) and `CANX KO` (`0.269038 < 1.111691` VX-661 & `1.041286` VX-445 `< 2.978676` Combo), and shuffled-predictor negative controls (`-0.003609` Arm A, `-0.126410` Arm B).
- **Arm A (`canx_chaperone_coupling`) — REFUTED:**
  - `baseline_dmso_ols` holdout MAE (`log2`): `0.451487`
  - `candidate` (`log2(Parental DMSO) + log2(CANX KO DMSO)`): `0.430771` (`+4.59%`, `relative_mae_reduction = 0.045884 < 0.10`)
  - 59-cluster bootstrap CI: 95% `[0.000723, 0.041823]`, **Bonferroni 97.5% `[-0.001739, 0.044782]`** (`ci975_low_gt_zero = false`). High collinearity (`r = 0.992301`) yields large opposing slopes (`+1.688711` vs. `-1.001083`) and does not reach the `10%` or Bonferroni threshold.
- **Arm B (`canx_ko_single_corrector_selectivity`) — REFUTED:**
  - `constant_dev_median` (`-0.108261`, locked in 5-fold dev CV `0.186444` vs. `canx_ko_dmso_ols` `0.191943`): holdout MAE `0.186365`
  - `candidate_partial_ridge` (`selected_alpha = 100.0`): holdout MAE `0.189235` (`-1.54%`, `relative_mae_reduction = -0.015397`, **Bonferroni 97.5% CI `[-0.038826, 0.036407]`**).

Both refuted hypotheses are logged in `failures.jsonl` and `ledger.json`.

## Running Verification

```bash
python -m unittest discover -s benchmarks/cftr_dms_canx_selectivity -p "test_*.py" -v
python benchmarks/cftr_dms_canx_selectivity/cftr_dms_canx_selectivity.py --check
python benchmarks/cftr_dms_canx_selectivity/independent_audit.py --check
```

When local git-ignored `.xlsx` files exist under `sources/downloads/tedman-2025/`, `--check` and `independent_audit.py --check` verify full raw `.xlsx` dual-sheet parsing and Cholesky $L L^T$ numerical parity (`tol <= 1e-6`). In CI without raw downloads, they verify `done_sha256`, `split_manifest_sha256`, codon disjointness, `0` `PRED-0010` holdout overlap, and verdict arithmetic.
