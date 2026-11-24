# Modulator-Free vs Trikafta-Bound $\Delta\text{F508}$ Cryo-EM Rescue Benchmark (`PRED-0014`)

> **Research benchmark only; not clinical or diagnostic guidance.** Zero patient or PHI data. Binary `*_responsive` fields record US FDA Structured Product Label (SPL) listing concordance (`listed_on_fda_spl` vs `not_listed_on_fda_spl`), not measured laboratory FRT nonresponse (`experimental_rescue_status = "unknown"`).

## What this benchmark tests

`PRED-0010` refuted standalone 1D AlphaMissense on TRIKAFTA v23 Table 6, and `PRED-0013` refuted single-state WT `5UAK` packing compliance. `PRED-0014` pre-registered (commit `9a23e16` in `done.md`) a two-state $\Delta\text{F508}$ Cryo-EM backbone conservation score:

$$s_{\text{stab}}(i) = \text{lDDT}_{\text{sym}}(i;\, \text{8EJ1},\, \text{8EIQ})$$

The score uses Mariani et al. (2013) published defaults ($R_0 = 15.0\text{ \AA}$, thresholds $\{0.5, 1.0, 2.0, 4.0\}\text{ \AA}$). It compares dephosphorylated modulator-free human $\Delta\text{F508}$ CFTR (PDB `8EJ1`, $6.9\text{ \AA}$, `1,127` Chain A $C_\alpha$) against phosphorylated Trikafta+ATP/Mg-bound human $\Delta\text{F508}$ CFTR (PDB `8EIQ`, $3.0\text{ \AA}$, `1,162` Chain A $C_\alpha$; Fiedorczuk & Chen, *Science* 2022, DOI `10.1126/science.ade2216`). Across $S_{\text{common}} = 1,121$ shared residues, it excludes covalent neighbors $|i - j| \le 1$ and the $\Delta\text{F508}$ `507-509` trans-peptide bond. It tests for $\ge 0.75$ ROC-AUC on an AST-blinded 24-variant untouched holdout panel (`holdout_open_annotations_2026-11-23.json`).

## Measured outcome (`refuted`)

Both independent scripts (`cryoem_f508del_rescue_benchmark.py` in `float64` and AST-isolated `independent_audit.py` in square-root-free `fractions.Fraction` arithmetic) agree to six decimal places:

| Metric on 24-variant panel (TRIKAFTA v23 Table 6) | Reduced Fraction | Decimal | Threshold / Baseline |
| --- | --- | --- | --- |
| Placed missense conditional ROC-AUC ($11/14$ missense in $S_{\text{common}}$) | `7/10` (`21/30` pairs) | `0.700000` | `< 0.75` (refuted) |
| All-14 missense abstention-aware ROC-AUC ($s_{\text{stab}}$) | `61/98` | `0.622449` | Beats standalone $s_{\text{AM}}$ (`22/49 = 0.448980`), `< 0.75` |
| All-24 class-gated abstention-aware ROC-AUC ($g_{\text{stab}}$) | `185/256` | `0.722656` | Ties $g_{\text{AM}}$ (`185/256 = 0.722656`), below class baseline (`95/128 = 0.742188`) |

Why `PRED-0014` was refuted:
1. **Unmodeled regions in $6.9\text{ \AA}$ apo $\Delta\text{F508}$ (`8EJ1`):** Three of 14 holdout missense variants abstain to `None` (`D110H`, `D648V`, and `E822K`). Under abstention-aware scoring (`None -> 0.5` tie), `D110H` and `E822K` tie Class I and canonical-splice negatives (`0.0`). That lowers all-24 $g_{\text{stab}}$ below the sequence-consequence class baseline (`0.722656` vs `0.742188`).
2. **Rigid Lasso/NBD1 scaffolding vs side-chain chemistry:** Among the 11 placed missense variants, unlisted Lasso mutants `W57G` (`0.877660`) and `D44G` (`0.847222`) sit in locally rigid segments. NBD1 `A613T` (`0.857143`, unlisted on TRIKAFTA v23 yet listed on ALYFTREK v9) also scores high and inverts against `H1085R`, `H139R`, and `S13F` (`9/30` pairs).

Following evaluation, `holdout_open_annotations_2026-11-23.json` carries `"panel_exposure_status": "exposed_development_reference"`. Future hypotheses require a fresh untouched holdout set.

## Offline verification

```powershell
python benchmarks/cryoem_f508del_rescue/cryoem_f508del_rescue_benchmark.py --check --falsify
python benchmarks/cryoem_f508del_rescue/independent_audit.py --check
python -m unittest discover -s benchmarks/cryoem_f508del_rescue -v
```
