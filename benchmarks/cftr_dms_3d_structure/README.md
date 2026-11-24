# `benchmarks/cftr_dms_3d_structure/` — Pure 3D Cryo-EM Structure-to-DMS Rescue Benchmark

Evaluates whether 3D structural descriptors derived offline from RCSB PDB Cryo-EM models (`6MSM`, `5UAK`, `8EIQ`; CC0 1.0) improve held-out prediction of CFTR modulator-rescued plasma-membrane expression (`log2(VX-661+VX-445 Mut/WT)`) beyond untreated baseline expression (`log2(DMSO Mut/WT)`) on the locked Tedman et al. (*eLife* 2025, `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`, CC BY 4.0) residue-clustered split.

> **Scope & biological boundaries:** Tedman et al. values measure HEK293T cell-surface plasma-membrane immunostaining (`Average Relative Intensity (Mut/WT)`). Surface abundance is **not** CFTR chloride conductance, patient lung function, or clinical guidance. `6MSM` `B_iso_or_equiv` is used strictly as a local Cryo-EM map-resolution/flexibility proxy, **not** as a physical atomic displacement parameter.

## Pre-Registered Protocol & Coordinator WT Parity Ruling

- **Pure 3D isolation:** Both `sub_domain` and `variant_class` indicators are excluded from baselines and candidate so `delta_MAE` reflects purely the 6 Cryo-EM features (`rbf_wjx`, `rbf_cv6`, `cn_10A_6msm`, `bfactor_6msm`, `disp_5uak_6msm`, `struct_unresolved`).
- **Dynamic `_atom_site.auth_seq_id` indexing:** Chain A C-alpha coordinates are keyed strictly by integer `_atom_site.auth_seq_id` (not `label_seq_id`, because `8EIQ` carries the `deltaF508` deletion and shifts `label_seq_id` by `-1` starting after codon `507`).
- **Coordinator WT amino-acid parity ruling (applied before committing `done.md`):** Across all `234` variants, `0` holdout variants and `0` unexpected `dev` variants fail WT parity against `6MSM`/`8EIQ` (aside from the `F508` deletion in `8EIQ`). Two `dev` source-table discrepancies were documented and handled cleanly:
  1. **Row `89` (`dev`, `M470V` / `1408A>G`, `codon=470`):** Tedman et al. 2025 uses the `M470` polymorphic background (`rs213950`), whereas canonical UniProt `P13569`, `6MSM`, `5UAK`, and `8EIQ` carry `V470` (`VAL`); `codon=470` remains mapped to `auth_seq_id=470`.
  2. **Row `217` (`dev`, `Q966=` / `2988G>A`, `codon=966`):** `Q966=` conflicts with `THR966` in `6MSM`, `5UAK`, `8EIQ`, and `P13569` (`2988G>A` corresponds to `Q996=`). Row `217` stays in `dev` (`234/234` row count preserved, zero speculative remapping) and is marked unmapped (`struct_unresolved = 1.0`, all `5` continuous 3D features imputed strictly inside training folds).

## Pre-Holdout Controls & Single-Shot Holdout Outcome (`results_2026-11-24.json`)

- **Pre-holdout controls (all passed):**
  - `TMD1` (`77..360`) union `TMD2` (`859..1155`) superposition: `8EIQ -> 6MSM` (`565` C-alpha) `RMSD_TMD = 0.843794` A (`< 1.5` A); `5UAK -> 6MSM` (`556` C-alpha) `RMSD_TMD = 5.009998` A; `WJX`–`CV6` min heavy-atom distance `30.013875` A (`36.894394` A centroid).
  - `#396` baseline parity on `dev` and `holdout` (`< 1e-6`): `constant_dev_median = 1.641988` `log2`; `dmso_ols`: `b0 = 2.21154`, `b1 = 0.865568`, `dev` 5-fold CV MAE `0.443911` `log2`.
  - Dev codon-permuted 3D null control (`seed = 20261124`): `0.445431` `log2` (`-0.3426%` vs `dmso_ols` CV MAE, `< 10%`).
- **Candidate `dev` CV vs. single-shot holdout (`67` variants, `59` codons):**
  - `dev` 5-fold CV (`selected_alpha = 100.0`): MAE `0.399628` `log2` (all `dev` VIFs `<= 1.38`).
  - `holdout` MAE (`log2` `VX-661+VX-445` `Mut/WT`): `constant_dev_median` **`1.077948`** (`1.593272` linear), `dmso_ols` **`0.451487`** (`0.541184` linear), `dmso_3d_cryoem_ridge` **`0.450918`** (`0.537163` linear).
- **Pre-registered decision gate (`refuted`):** Relative $\log_2$ MAE reduction over `dmso_ols` is **`+0.1261%`** (`0.001261 < 0.10`), with `10,000`-resample `59`-cluster 95% bootstrap CI **`[-0.065011, 0.066022]`**. Non-gated holdout breakdown shows minor gains in `NBD2` (`+0.109595`), `NBD1` (`+0.031672`), `R Domain` (`+0.033986`), and `Lasso` (`+0.397088`, `N=2`) offset by small losses in `TMD1` (`-0.029282`) and `TMD2` (`-0.065405`). Recorded honestly in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`.

## CI vs. Local Raw Data Mode

Tracked Cryo-EM mmCIF files (`tools/cftr-model/data/*.cif.gz`) enable full offline structural parity verification (`6MSM`/`5UAK`/`8EIQ` SHA-256 digests, Horn 4x4 quaternion superposition, `WJX`–`CV6` geometry, and WT sequence parity) even in CI when git-ignored `.xlsx` downloads are absent. When `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` exists locally under `sources/downloads/tedman-2025/`, `--check` and `independent_audit.py --check` also verify full raw `.xlsx` parsing and Cholesky $L L^T$ numerical parity (`<= 1e-6`).

## Run Checks

```bash
python -m unittest discover -s benchmarks/cftr_dms_3d_structure -p 'test_*.py' -v
python benchmarks/cftr_dms_3d_structure/cftr_dms_3d_structure.py --check
python benchmarks/cftr_dms_3d_structure/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/cftr_dms_3d_structure/ledger.json
```
