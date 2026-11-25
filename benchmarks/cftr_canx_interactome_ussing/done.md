# Done file: CFTR Calnexin Interactome & FRT `siCANX` Ussing Asymmetry Benchmark (`cftr_canx_interactome_ussing`)

Pre-registered and committed as a standalone commit on `feat/cftr-canx-interactome-ussing` before implementing evaluation code, running statistical tests, inspecting Arm 1 or Arm 2 outcomes, or writing `ledger.json`.

**Domain boundary:** This benchmark analyzes macroscopic Fischer Rat Thyroid (FRT) Ussing-chamber short-circuit current (`Delta Isc`) summaries and HEK293T TMT-18plex CFTR interactome proteomics (`log2` fold changes over the parental `C37` background) from Tedman et al. 2025 (*eLife*, DOI `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`, `CC BY 4.0`). Macroscopic epithelial `Delta Isc` and HEK293T co-immunoprecipitation fold changes are **not** single-channel patch-clamp burst durations and are **not** clinical human patient data or medical advice.

## Success test

The benchmark is complete when all deterministic checks exit `0` under the fixed protocol below without post-hoc retuning of inclusion rules, duplicate-handling rules, module definitions, or statistical thresholds:

1. **SHA-256-verified raw workbook provenance (git-ignored local files):**
   - Load workbooks solely from git-ignored `sources/downloads/tedman-2025/` (with clean relative worktree fallback; zero `.xlsx` copies committed to git) and verify SHA-256 digests on load:
     - `Copy of TMT_Labels.xlsx` (`5b16d08dfefce36878149259d823604d2153df8b027a48bd8cf69c3f1af5b41e`)
     - `Copy of CFTR_Variants_LogFoldChanges.xlsx` (`492df3ef266f27a5e18e1d9ce61b13356a8101b3e85d89c38ec87f5ca7208f7a`)
     - `Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx` (`d42ae782fe85163909396ea5ab4ee0915eec346e8ae420e419e80b5dd5e270c0`)
     - `Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx` (`a1ad7e325bab60ba6fe15e001fdddfd4fd14918d571e53ca45844c2f58f1ac1a`)
     - `Copy of FRT-WT_siCANX_Ussing-summary-for-paper.xlsx` (`bddaffb3b68449184b518a923c3c0d18ce51ed409662915d7b16fe7f8fedc628`)
   - Non-gated DMS concordance context optionally reads `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx` (`2a702ee3744b8039a9c7c4c74aec9bb30a91301b138d4948e0e2c81a1931e51b`, Row `213` `F508del`).

2. **Arm 1 (`frt_ussing_sicanx_asymmetry`) — FRT Ussing `Delta Isc` Protocol:**
   - **OOXML parsing:** Parse raw per-replicate `(Fsk, VX-770, Inh172)` cells from sheet `Raw-data-analyzed` of both Ussing workbooks via parsed cell-reference coordinates (`c.attrib["r"]`), ignoring author AVG/SEM/t-test summary blocks and derived sheets (`acute-treatments-analyzed-alone`, `percent-of-WT-analyses`). Read replicate columns `Cols C..L` in `FRT-F508del` (numbered `1..10` ahead of `M1 == 'AVG'`) and `Cols C..R` in `FRT-WT` (numbered `1..16` ahead of `S1 == 'AVG'`; all `10` non-`NT` groups occupy `Cols C..J`, within `Cols C..L`, while `NT` occupies `15` replicates across `C..O, Q..R`). Fail closed on any partial replicate or non-numeric replicate cell. Compute per-replicate `total_act = Fsk + VX-770`.
   - **Primary replicate rule (strict deduplication):** Within each sample group, detect and deduplicate replicate columns whose `(Fsk, VX-770, Inh172)` 3-tuple is bit-identical, keeping the earliest column and dropping subsequent duplicates (`Col C == Col D` in `7` of `10` non-`NT` `FRT-F508del` groups, reducing `N=4` to `N=3`; `FRT-WT` has zero duplicates, `N=6..10` non-`NT` and `N=15` `NT`). Record every detected duplicate pair neutrally and factually under `data_integrity_audit.duplicate_replicate_columns`.
   - **Primary gated hypotheses (Bonferroni `alpha = 0.025` per sub-test; no cross-arm OR-gate; `10,000`-resample stratified bootstrap with `seed = 20261125`):**
     - *Arm 1a (`f508del_combo_sicanx_drop`):* Compare `FRT-F508del siCANX, 661+445` vs. `FRT-F508del siNS, 661+445` on `total_act = Fsk + VX-770`:
       - Full two-sided permutation enumeration over all `C(n1 + n2, n1)` splits using `|mean(siCANX) - mean(siNS)|`: require `perm_p_two_sided <= 0.025`.
       - Stratified `10,000`-resample bootstrap (`Random(20261125)`) two-sided `97.5%` quantile CI (`[1.25%, 98.75%]`) on `ln(mean(siCANX) / mean(siNS))`: require `log_ratio_ci975_high < 0.0`.
       - Report Welch's `t` statistic and two-sided p-value as companion parametric context.
     - *Arm 1b (`sicanx_proportion_f508del_vs_wt`):* Compare proportional `siCANX` change in `FRT-F508del 661+445` vs. `FRT-WT DMSO` via difference-in-differences `did_log_ratio = ln(mean(F508del siCANX 661+445) / mean(F508del siNS 661+445)) - ln(mean(WT siCANX DMSO) / mean(WT siNS DMSO))`:
       - Full stratified difference-in-differences permutation enumeration over all `C(n_F_c + n_F_s, n_F_c) * C(n_W_c + n_W_s, n_W_c)` within-cell-line label swaps (`450,450` splits under deduplication): require one-sided `did_perm_p_one_sided <= 0.025` (`count(did_perm <= did_obs) / total_splits`; two-sided p-value reported as context).
       - Stratified `10,000`-resample bootstrap (independent `Random(20261125)` stream) one-sided `97.5%` upper quantile (`q_975` on `did_log_ratio < 0`): require `did_ci975_one_sided_high < 0.0` (two-sided `[1.25%, 98.75%]` CI reported as context).
     - **Arm 1 verdict rule:** Marked `supported` (`standing: "open"`) only if **both** Arm 1a and Arm 1b pass their `alpha = 0.025` permutation and `97.5%` bootstrap CI gates; otherwise marked `refuted` (`standing: "refuted"`).
   - **Non-gated sensitivity context:** Report identical Arm 1a and Arm 1b checks on as-deposited un-deduplicated columns (`sensitivity_as_deposited_n4`), and Spearman rank concordance across shared treatments (`DMSO`, `VX-661`, `VX-445`, `661+445`) vs. HEK293T DMS `F508del` (`Row 213` of `ALT 9-6-24 CFTR DMS Expanded Library Master Tables.xlsx`: 3 populated Parental treatments `DMSO`, `VX-445`, `661+445` because Parental `Col J` `VX-661` is blank `0/234`, plus 4 populated `All_variants CANX KO` treatments).

3. **Arm 2 (`canx_ko_proteostasis_remodeling`) — TMT-18plex Proteomics Protocol:**
   - **OOXML parsing:** Parse `Copy of CFTR_Variants_LogFoldChanges.xlsx` (`Sheet1`, `373` non-empty data rows `R03..R375`) and `Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx` (`Sheet1`, `236` non-empty data rows `R03..R238`), normalizing `\xa0 -> " "` in headers and failing closed if `> 5` non-numeric cells occur per workbook or if any non-numeric cell is other than `#DIV/0!`. Exclude normalization anchor row `CFTR` (`all 0.0`) before all statistical tests and apply pairwise-complete exclusion on any `#DIV/0!` cell (`RAB5A` in `I1366N (-CX)`; `LUZP1`, `POLE3` in `ΔF (-CX)`).
   - **Four canonical proteostasis modules (`.strip().lower()` on `Pathway`):**
     - `folding`: `{"folding", "protein folding"}`
     - `ups`: `{"ubiquitin proteasome system"}`
     - `autophagy_lysosome`: `{"autophagy-lysosome pathway"}`
     - `trafficking`: `{"trafficking"}`
   - **Primary gated hypotheses (2 separate BH-FDR `q < 0.05` families; no cross-family OR-gate):**
     - *Family 2a — Drug attenuation (`4` tests, BH-FDR `q < 0.05` across the 4 tests):* Paired two-sided Wilcoxon signed-rank test on `d_i = |LFC(DelF (-CX) + combo)| - |LFC(ΔF (-CX))|` per canonical module (excluding `d_i == 0`; combinatorial integer-doubled rank-sum dynamic programming when `n_nonzero <= 25`, tie-corrected normal approximation otherwise). Signed shift `LFC(combo) - LFC(DMSO)` reported as non-gated context.
     - *Family 2b — Variant enrichment (`20` tests, BH-FDR `q < 0.05` across the 20 tests):* Two-sided Mann-Whitney U test (tie-corrected normal approximation) comparing `LFC` of each canonical module against non-module background across all `5` variants (`WT(-CX)`, `V232D (-CX)`, `P5L (-CX)`, `T1036N (-CX)`, `I1366N (-CX)`). Unsigned `|LFC|` MWU tests reported as non-gated context.
     - **Arm 2 verdict rule:** Marked `supported` (`standing: "open"`) only if **both** `>= 1` canonical module achieves BH-FDR `q < 0.05` in Family 2a (`drug_attenuation`) **and** `>= 1` canonical module achieves BH-FDR `q < 0.05` in Family 2b (`variant_enrichment`); otherwise marked `refuted` (`standing: "refuted"`).
   - **Non-gated sensitivity context:** Evaluate `broad_proteostasis_sensitivity` defined as the union of the 4 canonical modules and `{"er proteostasis", "cytonuclear proteostasis", "mitochondrial proteostasis", "nuclear proteostasis", "extracellular proteostasis", "er", "er, trafficking", "vesicle trafficking", "quality control of membrane proteins"}`.

4. **Independent audit & deterministic repository gates:**
   - `independent_audit.py` (importing nothing from `cftr_canx_interactome_ussing.py` or `tools/`) validates `done_sha256`, independently parses all 5 `.xlsx` sheets when present locally, and checks all Arm 1 and Arm 2 computations within `tol <= 1e-6`.
   - Unit and planted-defect tests (`test_*.py`), `--check` (offline CI mode when `.xlsx` files are absent plus full parity locally), `independent_audit.py --check`, `check_claim_ledger.py --git`, `check_repo_organization.py`, `check_claim_strength.py`, `check_codeql_patterns.py`, and `check_anti_hallucination.py` all exit `0`.

## Positive control

Before interpreting Arm 1 or Arm 2 hypothesis outcomes, all technical and directional positive/negative controls must pass (aborting before outcome reporting if any fails):

1. **Arm 1 `Inh172` reversal parity:** Every group (`11` in `FRT-F508del` and `11` in `FRT-WT`) has `mean(Inh172) < 0` and `mean(Fsk + VX-770) > 0`, and across the `11` group means in each workbook `Pearson r(mean(Fsk + VX-770), -mean(Inh172)) > 0.95`.
2. **Arm 1 WT vs. F508del dynamic range:** `FRT-WT NT mean(Fsk + VX-770) > 5 * FRT-F508del NT mean(Fsk + VX-770)` (verified on both `Cols C..R` `N=15` and `Cols C..L` `N=10` for `FRT-WT NT`).
3. **Arm 1 F508del modulator rescue:** `FRT-F508del siNS, 661+445 mean(Fsk + VX-770) > FRT-F508del siNS, DMSO mean(Fsk + VX-770)`.
4. **Arm 2 `CFTR` normalization anchor, dispersion, & TMT design parity:** Row `CFTR` equals `0.0` (`|val| <= 1e-12`) across all `9` proteomics columns (`5` Variants + `4` Drug), every non-`CFTR` column has sample `std > 0.05`, and `Copy of TMT_Labels.xlsx` verifies `6` TMT-18plex runs (`18` labels `126..135N`).
5. **Arm 2 pathway-label permutation negative control (`10,000` draws, `seed = 20261125`):** Shuffling `Pathway` labels across non-`CFTR` proteins verifies calibrated per-family null false-discovery rate (`<= 0.08` at BH-FDR `q < 0.05` across Drug `4` tests and Variants `20` tests).

## Stopping rule

Single-shot evaluation on the 5 SHA-256-verified Tedman et al. 2025 workbooks after committing this `done.md`. Zero post-hoc modification of replicate deduplication rules, `Pathway` module sets, `#DIV/0!` handling, permutation seeds (`20261125`), or Bonferroni/BH-FDR thresholds after viewing Arm 1 or Arm 2 outcomes. Any control failure aborts before outcome reporting; any arm failing its pre-registered gate is recorded honestly as `refuted` in `results_2026-11-25.json`, `ledger.json`, `failures.jsonl`, and `README.md`.
