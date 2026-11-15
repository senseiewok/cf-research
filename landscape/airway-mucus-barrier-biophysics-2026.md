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

In `tools/discovery/mucus_barrier_screener.py`, the Stokes-Einstein-Ogston-Brenner (`SEOB`) steric and electrostatic hindrance kernel was benchmarked against 5 published CF sputum multiple-particle-tracking (`MPT`) carrier regimes (`EMPIRICAL_MPT_BENCHMARK_CASES`). Across the 5-regime empirical CF sputum MPT benchmark, all 5 published carrier regimes matched their empirical sputum diffusivity ratio intervals and expected feasibility verdicts (5 of 5 passed, 3 of 3 falsification controls passed):

| Regime ID | Carrier Description | $d_p$ (nm) | $\zeta$ (mV) | Mesh $\xi$ (nm) | Layer $L$ ($\mu\text{m}$) | Computed $D_{\text{eff}}/D_w$ | Empirical $D_{\text{eff}}/D_w$ Interval | Transit $\tau$ (s) | Deterministic Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `mpp_dense_peg_80nm` | Dense-PEG mucus-penetrating nanoparticle | 80 | -2.0 | 140 | 15 | 0.052478 | 0.025 to 0.12 | 724.84 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `mpp_hybrid_lipid_polymer_100nm` | Dense-PEG hybrid lipid-polymer LNP | 100 | -3.0 | 140 | 15 | 0.036923 | 0.015 to 0.08 | 824.87 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `aav_capsid_25nm_neutral` | Neutral AAV capsid vector | 25 | -3.5 | 140 | 20 | 0.175382 | 0.08 to 0.35 | 120.23 | `FEASIBLE_MUCUS_PENETRATING_CANDIDATE` |
| `uncoated_cationic_lipoplex_110nm` | Uncoated cationic lipoplex (0 percent PEG) | 110 | +28.0 | 140 | 20 | 0.000018 | 0.000001 to 0.0005 | $> 900$ | `REJECT_MUCOCILIARY_CLEARANCE_BEFORE_EPITHELIUM` |
| `oversized_conventional_particle_300nm` | Oversized conventional particle ($r_p = 150\text{ nm}$, `150 / 140 = 1.07` ratio) | 300 | -4.0 | 140 | 20 | 0.000516 | 0.00005 to 0.002 | $> 900$ | `REJECT_MUCOCILIARY_CLEARANCE_BEFORE_EPITHELIUM` |

---

## 3. Mechanistic Synthesis: Electrostatic Adhesion vs. Steric Sieving

Comparing the 100 nm dense-PEG hybrid nanoparticle (`D_eff/D_w` of 0.036923) against the 110 nm uncoated cationic lipoplex (`D_eff/D_w` of 0.000018), we infer that surface charge shielding dominates over modest diameter differences near 100 nm in CF sputum. Conversely, when particle radius $r_p = 150\text{ nm}$ ($d_p = 300\text{ nm}$) exceeds the untreated 140 nm sputum pore size (`150 / 140 = 1.07`), steric obstruction alone reduces `D_eff/D_w` to 0.000516 even with a near-neutral zeta potential (-4.0 mV), preventing epithelial arrival within the 900-second mucociliary clearance window.
