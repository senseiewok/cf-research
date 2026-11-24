# Done file: Pure 3D Cryo-EM Structure-to-DMS Rescue Benchmark (`cftr_dms_3d_structure`)

Pre-registered and committed as a standalone commit before running outcome evaluation, reading holdout intensities, or recording `ledger.json`.
Research benchmark only: plasma-membrane expression in HEK293T DMS assays is not channel conductance, clinical lung response, or medical advice.

## Success test

The benchmark is complete when all deterministic checks pass on the fixed protocol below without post-hoc retuning:

1. **Data provenance, locked split, and offline Cryo-EM coordinates:**
   - **DMS workbook:** Maintainer-ratified CC BY 4.0 file `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` (`sha256`: `2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b`, Tedman et al., *eLife* 2025, `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`) located in git-ignored `sources/downloads/tedman-2025/` (resolved across main checkout or `.worktrees/cf-research/...`, with graceful offline fallback in CI when absent). Parse sheet `All_variants` via sparse XML coordinates: `y_combo = log2(VX-661+VX-445 Mut/WT)` (Col `R`) and `x_dmso = log2(DMSO Mut/WT)` (Col `F`), requiring strictly positive finite values (`> 0`) and zero outcome imputation (`234/234` rows preserved).
   - **Locked residue-clustered split:** Reuse `benchmarks/cftr_dms_expression/split_manifest.json` (canonical `sha256`: `e6d60aa56d2ac87a9b82aaca8b3df5b560131943eb0ffb809f0834e8a9db3543`): `167` `dev` variants across `123` unique codons (including all `45` variants across the `30` exposed `PRED-0010` codons) and `67` `holdout` variants across `59` unique codons (`0` overlap with `dev` or `PRED-0010`).
   - **Tracked offline Cryo-EM mmCIFs (`tools/cftr-model/data/`, CC0 1.0):**
     - `6MSM.cif.gz` (`sha256`: `82933988abc3ed24a343819d410b39dc0d86913eaba4426dd945e3273e34fa59`) — WT phosphorylated ATP-bound CFTR
     - `5UAK.cif.gz` (`sha256`: `8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf`) — WT dephosphorylated ATP-free CFTR
     - `8EIQ.cif.gz` (`sha256`: `d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6`) — deltaF508 phosphorylated CFTR + Trikafta (`WJX`=VX-445, `CV6`=VX-661, `VX7`=VX-770)
   - **Dynamic mmCIF indexing:** Parse `_atom_site.*` columns dynamically and index Chain A C-alpha atoms (`group_PDB == "ATOM"`, `auth_asym_id == "A"`, `label_atom_id == "CA"`, `pdbx_PDB_model_num == "1"`) strictly by integer `_atom_site.auth_seq_id` (disallowing `_atom_site.label_seq_id`, because `8EIQ` carries the F508 deletion and shifts `label_seq_id` by `-1` starting after codon `507`). Verify codon `508` is `PHE` in `6MSM` and absent in `8EIQ`.

2. **Single primary arm & pure 3D structural isolation:**
   - Exclude `variant_class` and `sub_domain` indicators from both baselines and candidate so `delta_MAE` isolates purely 3D Cryo-EM structural features.
   - **Baselines (fitted on `dev`):** `constant_dev_median` (median of `y_combo` on `dev`) and `dmso_ols` (`y_combo = b0 + b1 * x_dmso`), locking `stronger_baseline` by lower 5-fold residue-clustered `dev`-CV MAE and verifying `< 1e-6` holdout MAE parity with `benchmarks/cftr_dms_expression/results_2026-11-23.json`.
   - **Candidate (`dmso_3d_cryoem_ridge`):** Partial-penalized Ridge regression (`0` L2 penalty on intercept `b0` and `x_dmso` slope `b1`; L2 penalty `alpha in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)` selected solely inside `dev` via 5-fold residue-clustered CV on `dev`-standardized 6-feature 3D vectors).

3. **Fixed 6-feature 3D vector & zero-leakage imputation:**
   - Superpose `8EIQ` and `5UAK` onto `6MSM` via Kabsch/SVD rigid alignment over common Chain A `auth_seq_id` C-alpha atoms in `TMD1` (`77..360`) union `TMD2` (`859..1155`). Apply the `8EIQ -> 6MSM` rigid transform to `8EIQ` `HETATM` heavy atoms (`type_symbol` not in `("H", "D")`) of `WJX` and `CV6`.
   - Features per codon `c`:
     1) `rbf_wjx = exp(-(d_wjx / 12.0)**2)` (`d_wjx` = minimum distance in Angstroms from `6MSM` C-alpha at `c` to superposed `WJX` heavy atoms)
     2) `rbf_cv6 = exp(-(d_cv6 / 12.0)**2)` (`d_cv6` = minimum distance in Angstroms from `6MSM` C-alpha at `c` to superposed `CV6` heavy atoms)
     3) `cn_10A_6msm` = count of other (`seq != c`) `6MSM` Chain A C-alpha atoms within `10.0` Angstroms
     4) `bfactor_6msm` = `6MSM` Chain A C-alpha `B_iso_or_equiv` (explicitly labeled as a local Cryo-EM map-resolution/flexibility proxy, not a physical atomic displacement parameter)
     5) `disp_5uak_6msm` = Euclidean distance (Angstroms) between TMD-superposed `5UAK` C-alpha and `6MSM` C-alpha at `c`
     6) `struct_unresolved` = `1.0` if `c` is absent from `6MSM` or `5UAK` Chain A `auth_seq_id` (or flagged unmapped per the Row `217` rule below), else `0.0`
   - Any missing continuous 3D feature is imputed strictly to the training-fold mean inside each `dev` CV fold, and to the full `dev` mean for `holdout` (zero cross-split or cross-fold leakage).

4. **One-shot holdout decision gate & independent verification:**
   - Evaluate `constant_dev_median`, `dmso_ols`, and `dmso_3d_cryoem_ridge` once on the locked `67`-variant (`59`-codon) `holdout`.
   - Compute a `10,000`-resample seeded (`seed=20261124`) residue-clustered bootstrap 95% confidence interval on `delta_MAE = MAE_stronger_baseline - MAE_candidate` over the `59` holdout codons.
   - **Pre-registered pass criterion (`supported`):** `relative_mae_reduction >= 0.10` and `bootstrap_ci95_low > 0.0`. Otherwise record `"refuted"` (or `"inconclusive"` if `relative_mae_reduction >= 0.10` while `bootstrap_ci95_low <= 0.0`) honestly in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`. Report non-gated per-sub-domain MAE breakdown and `dev` feature VIFs as context only.
   - Independent audit (`benchmarks/cftr_dms_3d_structure/independent_audit.py`, importing nothing from `cftr_dms_3d_structure.py` or `tools/`) independently solves rigid superposition via Horn's 4x4 quaternion Jacobi eigensolver and partial-penalized Ridge via Cholesky `L L^T`, checking `<= 1e-6` parity and all offline/CI invariants.
   - `unittest`, `--check`, `independent_audit.py --check`, `check_claim_ledger.py --git`, `check_anti_hallucination.py`, `check_repo_organization.py`, `tools/sources` tests, `check_codeql_patterns.py`, and `check_claim_strength.py` all exit `0`.

## Positive control

Before evaluating holdout outcomes, all four controls must pass (fail closed before touching holdout if any fails):

1. **WT amino-acid sequence parity (strict coordinator ruling):** Across all `234` variants, `label_comp_id` at `auth_seq_id` in `6MSM` and `8EIQ` must match the WT residue parsed from `Protein distinction`, with zero unexpected mismatches anywhere (`0` in `holdout`) and only these explicitly documented exceptions:
   - Codon `508` (`F508`) deletion in `8EIQ` (present as `PHE` in `6MSM` and `5UAK`, absent in `8EIQ`).
   - Row `89` (`dev`, `M470V` / `1408A>G`, `codon=470`): Tedman et al. 2025 carries the common CFTR `M470` polymorphic background (`rs213950`), whereas canonical UniProt `P13569`, `6MSM`, `5UAK`, and `8EIQ` carry `V470` (`VAL`); keep `codon=470` mapped to `auth_seq_id=470` for 3D feature extraction.
   - Row `217` (`dev`, `Q966=` / `2988G>A`, `codon=966`): Source-table discrepancy (`Q966=` vs `THR966` in `6MSM`, `5UAK`, `8EIQ`, and `P13569`, while `2988G>A` corresponds to `Q996=`); keep `split_manifest.json` (`codon_position=966`, `split=dev`) intact without exclusion or speculative remapping, and treat Row `217` as unmapped to a verified 3D residue (`struct_unresolved = 1.0`, all `5` continuous features set to missing and imputed strictly in-fold).
2. **TMD core RMSD & modulator pocket separation:** Assert `RMSD_TMD(8EIQ -> 6MSM) < 1.5` Angstroms over common Chain A `TMD1` union `TMD2` C-alpha atoms; record `8EIQ->6MSM` and `5UAK->6MSM` RMSDs plus `WJX`–`CV6` separation.
3. **`#396` baseline parity:** Assert `dmso_ols` `dev` fit and holdout MAE match `benchmarks/cftr_dms_expression/results_2026-11-23.json` within `1e-6`.
4. **Dev-only codon-permuted 3D null control (`seed=20261124`):** Permuting the 6-feature rows across unique `dev` codons yields `< 10%` 5-fold `dev`-CV MAE reduction relative to `stronger_baseline`.

## Stopping rule

Single-shot evaluation on the locked `67`-variant (`59`-codon) `holdout` after all pre-holdout controls pass. Zero retuning of features, RBF bandwidth (`12.0` Angstroms), neighbor cutoff (`10.0` Angstroms), TMD alignment ranges, imputation rules, Ridge `alpha` grid, or split boundaries after reading holdout outputs. Any control failure stops execution before reading holdout; any miss against `>= 10%` MAE reduction or `ci95_low > 0.0` is recorded as `refuted` (or `inconclusive`) in `results_2026-11-24.json`, `ledger.json`, and `failures.jsonl`.
