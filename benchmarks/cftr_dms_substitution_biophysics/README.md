# `benchmarks/cftr_dms_substitution_biophysics/` — Substitution Biophysics + UniProt `P13569` Topology DMS Rescue Benchmark

Evaluates whether 7 substitution-biophysics and sequence-topology descriptors computed from each variant's experimental WT amino acid and tracked UniProt `P13569` (`tools/cftr-model/data/p13569.json`, `1,480` aa, CC BY 4.0) improve held-out prediction of CFTR modulator-rescued plasma-membrane expression (`log2(VX-661+VX-445 Mut/WT)`) beyond untreated baseline expression (`log2(DMSO Mut/WT)`) on the locked Tedman et al. (*eLife* 2025, `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`, CC BY 4.0) residue-clustered split.

> **Scope & biological boundaries:** Tedman et al. values measure HEK293T cell-surface plasma-membrane immunostaining (`Average Relative Intensity (Mut/WT)`). Surface abundance is **not** CFTR chloride conductance, patient lung function, or clinical guidance.

## Pre-Registered Protocol & Variant Rules

- **Locked 7-feature vector (`dmso_subst_biophysics_ridge`, partial Ridge with `0` penalty on intercept & `x_dmso`):**
  1. `delta_blosum62`: `BLOSUM62(wt_exp, mut) - BLOSUM62(wt_exp, wt_exp)` (Henikoff & Henikoff 1992, *PNAS* 89:10915 `20x20` matrix)
  2. `delta_kd_hydropathy`: `KD(mut) - KD(wt_exp)` (Kyte & Doolittle 1982, *J. Mol. Biol.* 157:105)
  3. `abs_delta_charge`: `|charge(mut) - charge(wt_exp)|` (`K,R=+1`, `D,E=-1`, others `0`)
  4. `helix_breaker_in_helix`: `1.0` if `P13569` marks codon in `Helix` or `Transmembrane` and a missense component has `wt_exp in {'P','G'}` or `mut in {'P','G'}`, else `0.0`
  5. `tm_polarity_shock`: `1.0` if `P13569` marks codon in `Transmembrane` and a missense component has `|delta_kd| >= 3.0` or `abs_delta_charge >= 1.0`, else `0.0`
  6. `is_nonsense`: `1.0` for premature termination codons (`E60X`..`Q1313X`), else `0.0`
  7. `is_inframe_indel`: `1.0` for in-frame indels (`L138ins`, `F312DEL`, `I507del`, `F508del`), else `0.0`
- **Experimental-WT & multi-class handling (`234/234` variants):**
  - **Experimental-WT rule:** Substitution deltas use the variant's experimental WT amino acid (`M -> V` for `Row 89 M470V`), rather than canonical `P13569`'s `V470`.
  - **Compound missense (`Row 77 Q359K/T360K`):** Sums continuous substitution deltas across components (`delta_blosum62 = -10.0`, `delta_kd_hydropathy = -3.6`, `abs_delta_charge = 2.0`) and takes logical-OR across component positions for binary topology indicators.
  - **Synonymous (`E528E`, `T854T`, `Q966=`, `L1156=`):** All `7` features set to `0.0`.
  - **In-frame indels & nonsense:** Continuous substitution scales set to missing (`None`) and imputed strictly to the training-fold mean inside `dev` CV and to the full `dev` mean for `holdout` (`14` `dev` rows, `6` `holdout` rows).

## Pre-Holdout Controls & Single-Shot Holdout Outcome (`results_2026-11-24.json`)

- **Pre-holdout controls (all passed):**
  - `P13569` WT sequence parity across `234` variants (`235` components): `0` holdout mismatches, `0` unexpected mismatches, with the two documented `dev` exceptions (`row 89 M470V`, `row 217 Q966=`).
  - `#396` baseline parity (`<= 1e-6`): `constant_dev_median = 1.641988` `log2`; `dmso_ols`: `b0 = 2.21154`, `b1 = 0.865568`, `dev` 5-fold CV MAE `0.443911` `log2`.
  - Dev variant-permuted null biophysics control (`seed = 20261124`): `selected_alpha = 1000.0`, `permuted_dev_cv_mae_log2 = 0.446154` (`-0.5055%` vs `dmso_ols` CV MAE, `< 10%`).
- **Candidate `dev` CV, weights, and VIFs:**
  - `dev` 5-fold CV (`selected_alpha = 100.0`): MAE **`0.423982`** `log2` (`b0 = 2.197789`, `b1 = 0.854118`).
  - Fitted `dev` standardized weights: `delta_blosum62 = -0.028556`, `delta_kd_hydropathy = 0.061095`, `abs_delta_charge = -0.007885`, `helix_breaker_in_helix = 0.009452`, `tm_polarity_shock = 0.038402`, `is_nonsense = -0.119115`, `is_inframe_indel = -0.000721`.
  - `dev` VIFs: all below `10.0` (`max_vif_features_only = 1.193513`; `max_vif_with_x_dmso_log2 = 1.445349`; `0` flagged features).
- **Single-shot holdout (`67` variants, `59` codons) & decision gate (`refuted`):**
  - `holdout` MAE (`log2` `VX-661+VX-445` `Mut/WT`): `constant_dev_median` **`1.077948`** (`1.593272` linear), `dmso_ols` **`0.451487`** (`0.541184` linear), `dmso_subst_biophysics_ridge` **`0.414943`** (`0.541441` linear).
  - Relative `log2` MAE reduction over `dmso_ols`: **`+8.0942%`** (`0.080942 < 0.10`, `delta_MAE = 0.036544`), with `10,000`-resample `59`-cluster 95% bootstrap CI **`[-0.001247, 0.080596]`**. Recorded honestly as `refuted` (`REFUTED_BELOW_10PCT_MAE_REDUCTION_THRESHOLD`) in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`.
  - **Non-gated context comparison:** `+8.0942%` here vs `#396` domain/class Ridge (`+3.7001%`) and `#401` pure 3D Cryo-EM Ridge (`+0.1261%`). Per-`sub_domain` holdout breakdown shows the largest reduction in `TMD2` (`0.519267 -> 0.401503`, `delta = +0.117764`) and `TMD1` (`+0.018911`), with slight regressions in `NBD1` (`-0.034494`), `R Domain` (`-0.020109`), and `Lasso` (`-0.001846`, `N=2`).

## CI vs. Local Raw Data Mode

Tracked UniProt `P13569` data (`tools/cftr-model/data/p13569.json`) enables offline verification (`p13569.json` SHA-256 against `ATTRIBUTION.md`, BLOSUM62 symmetry, Kyte-Doolittle anchors, WT sequence parity across all `235` components, and `7`-feature extraction across all `234` manifest rows) even in CI when git-ignored `.xlsx` downloads are absent. When `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` exists locally under `sources/downloads/tedman-2025/`, `--check` and `independent_audit.py --check` also check raw `.xlsx` parsing and Cholesky `L L^T` numerical parity (`<= 1e-6`).

## Run Checks

```bash
python -m unittest discover -s benchmarks/cftr_dms_substitution_biophysics -p 'test_*.py' -v
python benchmarks/cftr_dms_substitution_biophysics/cftr_dms_substitution_biophysics.py --check
python benchmarks/cftr_dms_substitution_biophysics/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/cftr_dms_substitution_biophysics/ledger.json
```
