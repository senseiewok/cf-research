# `benchmarks/cftr_dms_expression/` — Quantitative CFTR Expression Benchmark (`Council Item #1`)

Empirical CFTR plasma-membrane expression benchmark built on the maintainer-ratified CC BY 4.0 dataset from **Tedman et al., *eLife* 2025** (`10.7554/eLife.107180`; Mendeley Data `10.17632/ksfmxv8zts.1`).

> **Scope & biological boundary:** Values measure HEK293T cell-surface CFTR plasma-membrane immunostaining (`Average Relative Intensity (Mut/WT)`). Surface expression is **not** channel Cl- conductance, patient lung function, or medical advice.

## Pre-Registered Findings (`results_2026-11-23.json`)

- **Cross-study positive control (`AGM ALT Library Comparison.xlsx`, `N = 103`):** Pearson $r = 0.998856$ (linear), $r = 0.998275$ ($\log_2$), Spearman $\rho = 0.965767$ between McKee et al. 2023 (`AGM`) and Tedman et al. 2025 (`ALT`).
- **Zero-leakage residue-clustered split:** All `234` parental variants (`All_variants`, zero missing exclusions on `DMSO` / `VX-661+VX-445`) were grouped strictly by codon; all `45` variants matching the `30` exposed `PRED-0010` codons were forced into `dev` (`167` variants, `123` codons), leaving `67` unexposed variants across `59` codons in `holdout`.
- **Holdout MAE comparison ($\log_2$ `VX-661+VX-445` `Mut/WT`):**
  - Constant `dev` median (`1.641988` $\log_2$): **`1.077948`** (`1.593272` linear).
  - Untreated expression OLS (`y_combo = 2.21154 + 0.865568 * x_dmso`): **`0.451487`** (`0.541184` linear).
  - Partial-penalized Ridge (`dmso` unpenalized, `alpha = 100.0` selected via 5-fold residue-clustered `dev` CV): **`0.434782`** (`0.573246` linear).
- **Pre-registered decision gate (`refuted`):** Adding CFTR sub-domain and variant-class indicators gave a `+3.70%` relative $\log_2$ MAE reduction over untreated-expression OLS, with a `59`-cluster 95% bootstrap CI of `[-0.028005, 0.058719]`, falling short of the pre-registered `10%` threshold and zero-excluding interval. Untreated baseline expression accounts for most of the variance on the multiplicative $\log_2$ scale.

## CI vs. Local Raw Data Mode

Raw downloaded `.xlsx` files remain git-ignored under `sources/downloads/tedman-2025/`. Locally `--check` and `independent_audit.py --check` re-parse the spreadsheets from disk and verify numerical parity; when git-ignored downloads are absent in CI, unit tests and split/hash checks still validate the protocol and sparse XML parser.

## Run Checks

```bash
python -m unittest discover -s benchmarks/cftr_dms_expression -p 'test_*.py' -v
python benchmarks/cftr_dms_expression/cftr_dms_benchmark.py --check
python benchmarks/cftr_dms_expression/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/cftr_dms_expression/ledger.json
```
