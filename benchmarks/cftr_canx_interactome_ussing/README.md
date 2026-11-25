# CFTR Calnexin Interactome Remodeling & FRT `siCANX` Ussing Asymmetry Benchmark

This directory implements a pre-registered two-arm empirical benchmark on the untapped **Tedman et al. 2025** (*eLife*, DOI `10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`, `CC BY 4.0`) proteomics and electrophysiology summary workbooks.

**Epistemic & clinical scope boundary:** Measurements reflect macroscopic Fischer Rat Thyroid (FRT) Ussing-chamber short-circuit current (`Delta Isc`) summaries and HEK293T TMT-18plex CFTR co-immunoprecipitation `log2` fold changes over parental `C37`. Macroscopic epithelial `Delta Isc` and HEK293T interactome fold changes are **not** single-channel patch-clamp burst durations, clinical human patient outcomes, or medical advice.

## Pre-Registered Protocol & Data Integrity Audit

`done.md` was committed prior to outcome evaluation (`done_sha256: d6eeee7cf61baf2cee4576f6094c278a623b214d3667effe42ff05c14360ca82`), reading the five SHA-256-verified `.xlsx` workbooks solely from git-ignored `sources/downloads/tedman-2025/`:

- **Ussing replicate deduplication audit (`data_integrity_audit`):** Parsing `Raw-data-analyzed` via parsed OOXML cell coordinates (`Cols C..L` in `FRT-F508del`, `Cols C..R` in `FRT-WT`, ignoring author summary blocks) detected **7 bit-identical `Col C == Col D` duplicate replicate 3-tuples** among the `10` non-`NT` `FRT-F508del` sample groups (`siNS, DMSO`; `siNS, VX-661`; `siNS, VX-445`; `siCANX, VX-661`; `siCANX, VX-445`; `siCANX, 661+445`; `siCANX, ETI`), reducing each from `N=4` to `N=3`. `FRT-WT` contained `0` duplicate replicate columns (`N=15` for `NT`, `N=6..8` for non-`NT`).
- **Proteomics quality audit:** Both `CFTR_Variants_LogFoldChanges.xlsx` (`373` rows) and `CFTR_F508del_Drug_LogFoldChanges.xlsx` (`236` rows) carry `CFTR = 0.0` as the normalization anchor (excluded from all module tests) and `3` `#DIV/0!` cells (`RAB5A` in `I1366N (-CX)`; `LUZP1`, `POLE3` in `ΔF (-CX)`), handled via pairwise-complete exclusion. `Copy of TMT_Labels.xlsx` matches a `6`-run `18`-plex (`108` channel entries) experimental design.

## Observed Results (Single-Shot Evaluation)

- **Arm 1 Positive Controls (passed):** `22/22` groups show `mean(Fsk + VX-770) > 0` and `mean(Inh172) < 0`; pooled `22`-group `Pearson r = 0.991056` (`FRT-F508del` `r = 0.998260`, `FRT-WT` `9` non-`ETI` groups `r = 0.993663`, all-`11` `FRT-WT` `r = 0.884189` due to chronic `ETI` pre-treatment basal current); `FRT-WT NT` (`79.227207` `N=15`; `85.022870` `N=10` `Cols C..L`) exceeds `FRT-F508del NT` (`1.733125`) by `45.71x` (`49.06x`); `FRT-F508del siNS, 661+445` (`35.393250`) exceeds `siNS, DMSO` (`1.737333`) by `20.37x`.
- **Arm 1 (`frt_ussing_sicanx_asymmetry`) — REFUTED (`alpha = 0.025` Bonferroni per sub-test):**
  - *Arm 1a (`f508del_combo_sicanx_drop`):* `FRT-F508del 661+445` `siCANX` (`24.226667`, `N=3`) drops `-31.55%` vs. `siNS` (`35.393250`, `N=4`), `ln` ratio `-0.379067` (`10,000`-resample bootstrap `97.5%` CI `[-0.632096, -0.100299]`), **but** full permutation enumeration (`3/35` splits) yields `p = 0.085714 > 0.025` (Welch `t = -2.746379`, `df = 3.583194`, `p = 0.058231`) -> **FAIL**.
  - *Arm 1b (`sicanx_proportion_f508del_vs_wt`):* `FRT-WT DMSO` also drops `-27.74%` under `siCANX` (`63.695250` vs. `88.149625`, `ln` ratio `-0.324926`), yielding `did_log_ratio = -0.054142`, one-sided stratified DiD permutation `p = 0.408714 > 0.025` (`450,450` splits; two-sided `0.813573`), and bootstrap one-sided `97.5%` upper bound `+0.224171 > 0` -> **FAIL**.
  - *Sensitivity context:* Retaining un-deduplicated `N=4` columns yields `FRT-F508del 661+445` `siCANX = 25.895000` (`-26.84%`, `perm p = 0.085714`, Welch `p = 0.044307`) and reverses `did_log_ratio` sign to `+0.012454` (`p_1s = 0.524227`). Across shared treatments (`DMSO`, `VX-661`, `VX-445`, `661+445`), FRT Ussing means match HEK293T DMS `F508del` rank order (`Spearman rho = 1.0` across `3` populated Parental DMS treatments and `4` `CANX KO` DMS treatments).
- **Arm 2 Controls (passed):** `CFTR = 0.0` across all `9` proteomics columns (`std = 0.322243..0.703448`); `10,000`-draw `Pathway` permutation null yields calibrated FWER `0.0578` (Drug) and `0.0432` (Variants) (`<= 0.08`).
- **Arm 2 (`canx_ko_proteostasis_remodeling`) — REFUTED (`BH-FDR q < 0.05` in both families):**
  - *Family 2a (`drug_attenuation`, `4` combinatorial DP Wilcoxon tests):* `0/4` pass `q < 0.05` (`folding` `n=10`, median `|LFC(combo)| - |LFC(DMSO)| = -0.163479`, `p = 0.160156`, `q = 0.640625`; `ups` `n=20`, `q = 0.701181`; `autophagy_lysosome` `n=8`, `q = 0.701181`; `trafficking` `n=8`, `q = 0.701181`).
  - *Family 2b (`variant_enrichment`, `20` MWU tests):* `0/20` pass `q < 0.05` (minimum `q = 0.694575` for `WT(-CX):folding` `p = 0.048141` and `I1366N (-CX):trafficking` `p = 0.070789`).
  - *Sensitivity context (`broad_proteostasis_sensitivity`):* `59` paired Drug proteins give `p = 0.691907`, and `107` broad proteostasis proteins across Variants give minimum `q = 0.499193`.

## Running Verification

```bash
python -B -m unittest discover -s benchmarks/cftr_canx_interactome_ussing -p "test_*.py" -v
python benchmarks/cftr_canx_interactome_ussing/cftr_canx_interactome_ussing.py --check
python benchmarks/cftr_canx_interactome_ussing/independent_audit.py --check
```

When local git-ignored `.xlsx` files exist under `sources/downloads/tedman-2025/`, `--check` and `independent_audit.py --check` re-parse all five workbooks independently (`xml.dom.minidom`, `Decimal`, bitmask/two-pointer permutations, `Counter` subset-sum Wilcoxon, and direct pairwise Mann-Whitney U) to `tol <= 1e-6`. In CI without `.xlsx` downloads, they verify `done_sha256`, `ledger.json`, BH-FDR step-up values, and decision-gate arithmetic.
