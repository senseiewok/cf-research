# Cross-Disease Mechanistic Bridge & Council Falsification Dossier (`L-0039`–`L-0042`, `T-0179`)

**Date:** 2026-10-28  
**Council Review:** Orthogonal Frontier Challenger + Local Second-Family Challenger (`gemma4:31b-it-q4_K_M`) + Deterministic Verification Suite (`check_leads.py`, `check_anti_hallucination.py`, `test_cross_domain_miner.py`)  
**Epistemic Tier:** Tier T2 (biophysical models with explicit `BiophysicalParameter` provenance and negative-control falsification gates; research leads registered in `leads/leads.json` as `L-0039` to `L-0042`)  
**Research Boundary:** Independent open-source computational research; not clinical guidance or treatment recommendations. Zero clinical efficacy extrapolation across diseases is permitted.

---

## 1. Plain-Language Overview ("Just If It Makes Sense")

Academic and clinical teams studying Cystic Fibrosis (`CF`) rarely have time to read across retina, liver, muscular dystrophy, or sister airway registries, yet several biological bottlenecks in CF obey the same physical and molecular laws studied in those sister fields. To make sure cross-disease exploration never dilutes the CF curative mission or overclaims clinical benefit, the Multi-Family Council established two hard invariants before any new cross-disease lead is admitted:

1. **Bidirectional CF Return (`INV-CF-ANCHOR`)**: Every cross-disease lead or model must start from and return directly to a specific CF curative bottleneck—specifically the $8.0\%$ (`CFFPR 2024`), $7.6\%$ (`CCFR 2024`), and $5.8\%$ adult / $11.5\%$ pediatric (`ACFDR 2025`) modulator-ineligible Class I/nonsense/splicing cohort, refractory missense domain-interface misfolding, dehydrated muco-obstructive mucus, or multidrug-resistant airway pathogen persistence.
2. **Mechanistic Isomorphism & Negative-Control Falsification (`INV-MECH-ISOMORPHISM`)**: A connection is valid only when governed by a shared mathematical equation or structural target class (`tools/discovery/cross_domain_miner.py` — `CrossDiseaseMechanisticBridge`), paired with an explicit negative control that rejects superficial clinical analogy.

---

## 2. Four Council-Verified Cross-Disease Bridges (`CDB-001` to `CDB-004`)

| Bridge ID & Lead | CF Curative Bottleneck (`INV-CF-ANCHOR`) | Sister Conditions / Targets | Shared Governing Law (`INV-MECH-ISOMORPHISM`) | Built-In Negative Control (Falsification Gate) |
| :--- | :--- | :--- | :--- | :--- |
| **`CDB-001` (`L-0039`)** | Airway mucociliary clearance & $15\text{ Hz}$ thixotropic liquefaction across all CF genotypes | Primary Ciliary Dyskinesia (`PCD`), Non-CF Bronchiectasis (`NCFB`), Alpha-1 Antitrypsin Deficiency (`AATD`) | Osmotic periciliary brush compression, viscoelastic yield stress $\tau_y(\phi_{\text{solids}}, C_{\text{DNA}})$, and registered `FEV1` / `LCI2.5` / exacerbation endpoint facets (`10.1056/NEJMra1813799`, `10.1038/s41598-025-16178-z`) | Recombinant human DNase I (`10.1073/pnas.87.23.9188`) lowers viscosity only when extracellular neutrophil DNA ($C_{\text{DNA}} \ge 1.5\text{ mg/mL}$) dominates purulent CF plugs; it is rejected by the model in low-DNA mucin-only plugs (`NCFB`). |
| **`CDB-002` (`L-0040`)** | Refractory Class II and rare missense `CFTR` (`ABCC7`) domain-interface misfolding & ERAD degradation | `ABCA4` (Stargardt retinopathy), `ABCB11` (`PFIC2`), `ABCC6` (`PXE`) | Boltzmann domain-interface folding stability $\Delta\Delta G_{\text{fold}}(\text{TMD–NBD})$ and `RNF5`/`CHIP` ubiquitin-ligase endoplasmic reticulum quality-control escape (`10.1126/science.ade2216`, `10.3389/fphar.2012.00201`, `10.1091/mbc.e14-04-0935`) | `CFTR` (`ABCC7`) is the sole ATP-gated ion channel in the human ABC superfamily; equating CFTR channel open-probability potentiation ($P_o$) with active ABC pump substrate turnover ($k_{\text{cat}}$) is deterministically rejected. |
| **`CDB-003` (`L-0041`)** | Class I nonsense `CFTR` alleles (`G542X`, `W1282X`, `R553X`, `R1162X`) in the modulator-ineligible cohort | Nonsense alleles in `DMD`, `USH2A`, and `IDUA` sharing `UGA`/`UAG`/`UAA` +4 nucleotide contexts | Coupled product law: $I_{\text{salvage}} = f_{\text{NMD\_escape}} \cdot P_{\text{readthrough}} \cdot P_{\text{fold}}(\Delta\Delta G_{\text{cavity}} + \Delta\Delta G_{\text{corrector}})$ (`10.1172/jci154571`, `10.1021/acs.biochem.5c00590`, `10.1016/j.jcf.2025.09.006`, `10.1038/s41586-023-06133-1`) | Standalone readthrough of a structural repeat (`DMD`) does not transfer to `CFTR` `W1282X`/`G542X` ($\Delta\Delta G_{\text{cavity}} \approx +2.8\text{ kcal/mol}$) without folding correctors (`trimodal_salvage_required = True`). |
| **`CDB-004` (`L-0042`)** | Full-length `CFTR` mRNA ($4{,}443\text{ bp}$) co-translational domain folding (`INV-tAI-PROFILE`) & dual-AAV split-intein assembly | `ABCA4` dual-AAV intein trans-splicing, Factor VIII (`F8`), and Dysferlin (`DYSF`) large-cassette reconstitution | Inter-domain ribosomal pause-valley correlation ($r \ge 0.85$) between `NBD1`, `R`, and `TMD2` (`10.1016/j.cell.2015.02.029`, `10.1126/scitranslmed.aav4523`, `10.1016/j.ymthe.2017.05.013`) | Naive global Codon Adaptation Index maximization ($\text{CAI} = 1.0$) erases inter-domain pause valleys and triggers co-translational domain collision (`COLLISION_RISK`). |

---

## 3. Reproducible Verification Commands

```powershell
python tools/leads/render_leads.py --check
python tools/leads/check_leads.py
python tools/discovery/cross_domain_miner.py --bridges
python -m unittest tools/discovery/test_cross_domain_miner.py -v
python tools/claims/check_anti_hallucination.py
```
