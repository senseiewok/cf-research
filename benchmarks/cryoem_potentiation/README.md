# Cryo-EM 3D Potentiation-vs-Folding Pocket Proximity Benchmark (`PRED-0011`)

Research benchmark only; not clinical, diagnostic, or prescribing guidance. Zero patient or PHI data.

## Question

Following `PRED-0010` (where 1D AlphaMissense sequence pathogenicity did not separate CFTR modulator-responsive variants on the 20-variant `NBD1`/`TMD2`/`NBD2_C_terminal` held-out split), does $E(3)$-invariant 3D Cryo-EM proximity to the ATP/ivacaftor potentiation machinery versus corrector/ICL4-folding pockets in `PDB 8EIQ` (`tools/cftr-model/data/8EIQ.cif.gz`, CC0 1.0) predict US FDA `KALYDECO` v36 Table 3 (`openfda-label`, NDA 203188) ivacaftor monotherapy responsiveness on the same fixed `2026-11-21.v2` held-out `test` split?

## Outcome (`PRED-0011`: `supported` on held-out `test`, with explicit `Lasso`/`TMD1` boundary)

- **Held-out `test` missense variants (`11` placed variants: `5` `KALYDECO`-responsive, `6` non-responsive)**:
  - `s_3D(i) = d_fold(i) / (d_pot(i) + d_fold(i))`: conditional ROC-AUC **`0.9333` (`14/15` = `28/30` pairs)** vs standalone `s_AM(i) = 1 - AM(i)` at `0.6667` (`2/3`) and constant-missense sequence baseline `0.5000`.
  - Only `2/30` (`1/15`) missense pairs were discordant: `TM12` pore-lining `D1152H` (`KALYDECO=1`, `s_3D = 0.4667`) scored below `NBD1` proline helix-breaker `L467P` (`KALYDECO=0`, `s_3D = 0.6085`) and `NBD2` `N1303K` (`KALYDECO=0`, `s_3D = 0.5689`).
- **All `20` held-out `test` variants (`7` `KALYDECO`-responsive, `13` non-responsive)**:
  - Class-gated `g_3D`: abstention-aware ROC-AUC **`0.8242` (`75/91` pairs)** vs deterministic sequence-consequence class baseline `0.6813` (`62/91`) and identically gated `g_AM` `0.7363` (`67/91`).
- **Documented biophysical scope boundary (`Lasso_N_terminal` / `TMD1` `train` split)**:
  - On the `7` `train` missense variants in `Lasso`/`TMD1`, static `s_3D` gave conditional ROC-AUC `0.2500` (`3/12` pairs) because mild residual-trafficking/conductance alleles (`P67L`, `R117H`, `L206W`) sit near transmembrane corrector cavities yet retain apical expression, while `TM6` `R347P` (`12.08 Å` from `VX7`) disrupts a pore salt bridge and fails ivacaftor monotherapy.
  - Per `docs/lab-focus-2026-11-14.md`, an external reader outside the lab still precedes calling any claim a result on a public page.

## Verification

```powershell
python -m unittest discover -s benchmarks/cryoem_potentiation -v
python benchmarks/cryoem_potentiation/cryoem_potentiation_benchmark.py --check --falsify
python benchmarks/cryoem_potentiation/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/cryoem_potentiation/ledger.json
```
