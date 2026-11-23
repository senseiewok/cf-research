# Airway Mucus Barrier Biophysics & Nanoparticle Transit Benchmark (2026)

*Track 2 Dossier — Airway Delivery & Mucus Barrier Biophysics (`T-0941` / `PRED-0004`)*
*Primary Sources:* `conte-2022-acsami`, `comegna-2021-scirep`, `alton-2024-erj-rsiv-fhn-lentivirus-preclinical`, `clinicaltrials-gov`
*Verification Sidecars:* `airway-mucus-barrier-biophysics-2026-claims.json`, `airway-mucus-barrier-biophysics-2026-evidence.txt`, `airway-mucus-barrier-biophysics-2026-evidence-map.json`, `airway-mucus-barrier-biophysics-2026-allow.txt`

---

## 1. Biophysical Constraints of the Cystic Fibrosis Airway Mucus Barrier

Every inhaled genetic or nucleic-acid therapy for Cystic Fibrosis—whether lipid nanoparticle (`LNP`) mRNA, base-editing ribonucleoprotein, adeno-associated virus (`AAV`), or pseudotyped lentivirus—must cross the dehydrated, hyper-viscoelastic airway mucus blanket before mucociliary clearance sweeps the periciliary airway surface (reference clearance window: 900 seconds or 15 minutes):

- **Hybrid Lipid-Polymer & Dense-PEG Shielding (`conte-2022-acsami`, DOI: `10.1021/acsami.1c14975`):** Conte G et al. 2022 reported hybrid lipid-polymer nanoparticles for pulmonary nucleic acid delivery with hydrodynamic diameters of 80 to 120 nm and dense polyethylene glycol (`PEG`) surface shells. Dense polyethylene glycol shells with near-neutral zeta potential between -5 mV and +2 mV shield nanoparticles from electrostatic adhesion to negatively charged mucin fibers in cystic fibrosis sputum. In vivo mucociliary transit clearance half-life values in human patients for the Conte G et al. 2022 formulation are `"not stated"` in the in vitro study record.
- **Epithelial Base Editing Delivery Requirement (`comegna-2021-scirep`, DOI: `10.1038/s41598-021-85549-z`):** Comegna M et al. 2021 demonstrated base editing restoration of `CFTR` function in human cystic fibrosis airway epithelial models, highlighting that non-viral carriers must traverse the viscoelastic mucus layer before reaching apical epithelial cells.
- **Viral Apical Entry Vector (`alton-2024-erj-rsiv-fhn-lentivirus-preclinical`, DOI: `10.1183/13993003.01683-2023`):** Alton EWFW et al. 2024 evaluated F/HN-pseudotyped `rSIV.F/HN` lentiviral vectors for cystic fibrosis airway gene transfer, demonstrating apical transduction of ciliated and basal airway epithelial cells in preclinical models.
- **Mucolytic Mesh Expansion (Dornase Alfa Pretreatment):** Untreated purulent cystic fibrosis sputum has an average mucin and extracellular DNA mesh pore size of 140 nm, whereas recombinant human DNase (dornase alfa) pretreatment expands the effective mesh pore size to 220 nm (`220 - 140 = 80` nm increase in effective mesh pore diameter).

---

## 2. 5-Regime CF Sputum Multiple-Particle-Tracking (MPT) Benchmark (`T-0941` / `PRED-0004`)

In `tools/discovery/mucus_barrier_screener.py`, the Renkin-Flory-Stokes-Einstein steric (`f_steric = (1 - d_p / d_pore)^2` for `d_p < d_pore`, else `0`) and electrostatic screening model outputs simulated quantities (`MODEL_VALUE`) benchmarked against 5 literature-defined CF sputum multiple-particle-tracking (`MPT`) carrier regimes (`EMPIRICAL_MPT_BENCHMARK_CASES`). All 5 simulated regimes fall within the target diffusivity-ratio bounds and match their expected classification verdicts (`5` of `5` passed, `3` of `3` falsification controls passed); these simulated screening values are not clinical or in vivo measurements and exclude macrophage clearance, sputum plug heterogeneity, airway impaction, endosomal escape, and clinical efficacy (`not_modeled_in_vivo`):

| Regime ID | Carrier Description | $d_p$ (nm) | $\zeta$ (mV) | Mesh $\xi$ (nm) | Layer $L$ ($\mu\text{m}$) | Simulated $D_{\text{eff}}/D_w$ (`MODEL_VALUE`) | Target Bounds | Simulated Transit $\tau$ (s) (`MODEL_VALUE`) | Deterministic Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `mpp_dense_peg_80nm` | Dense-PEG mucus-penetrating nanoparticle | 80 | -2.5 | 140 | 25 | 0.052478 | 0.02 to 0.15 | 724.84 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `mpp_hybrid_lipid_polymer_100nm` | Dense-PEG hybrid lipid-polymer LNP | 100 | -3.0 | 140 | 15 | 0.023324 | 0.01 to 0.08 | 733.90 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `aav_capsid_25nm_neutral` | Neutral AAV capsid vector | 25 | -4.0 | 140 | 25 | 0.192784 | 0.15 to 0.30 | 61.66 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `uncoated_cationic_lipoplex_110nm` | Uncoated cationic lipoplex (0 percent PEG) | 110 | +30.0 | 140 | 25 | 0.000054 | 0.0 to 0.005 | 975498.44 | `ELECTROSTATICALLY_TRAPPED_ON_MUCIN_FIBERS` |
| `oversized_conventional_particle_300nm` | Oversized conventional particle (`300 / 140 = 2.14` diameter-to-pore ratio) | 300 | -3.0 | 140 | 25 | 0.0 | 0.0 to 0.0 | `None` (sterically blocked) | `STERICALLY_TRAPPED_BY_SPUTUM_MESH` |

---

## 3. Mechanistic Synthesis: Electrostatic Adhesion vs. Steric Sieving

Comparing the simulated 100 nm dense-PEG hybrid nanoparticle (`D_eff/D_w` of 0.023324) against the simulated 110 nm uncoated cationic lipoplex (`D_eff/D_w` of 0.000054), we infer under this screening model that surface charge shielding produces a larger simulated diffusivity difference than the diameter difference between 100 nm and 110 nm in 140 nm sputum. Conversely, when simulated particle diameter $d_p = 300\text{ nm}$ meets or exceeds the modeled 140 nm sputum pore diameter (`300 / 140 = 2.14`), the screener's hard Renkin cutoff sets `f_steric = 0` and simulated `D_eff/D_w` to `0.0` even at a near-neutral zeta potential (-3.0 mV); this simplified single-pore threshold is a screening assumption and does not preclude residual transport through wider pores in heterogeneous human sputum.
