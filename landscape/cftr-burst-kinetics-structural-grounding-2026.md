# CFTR Single-Channel Burst Kinetics & Cryo-EM Structural Grounding (2026)

*Track 4 Dossier — CFTR Single-Channel Biophysics & Cryo-EM Structural Biology*
*Primary Sources:* `rcsb-pdb-6msm`, `rcsb-pdb-5uak`, `csanady-2006-jgp-irreversible-cycle`, `csanady-2010-pnas-burst-coupling`, `csanady-2019-physrev-cftr-review`
*Verification Sidecars:* `cftr-burst-kinetics-structural-grounding-2026-claims.json`, `cftr-burst-kinetics-structural-grounding-2026-evidence.txt`, `cftr-burst-kinetics-structural-grounding-2026-evidence-map.json`, `cftr-burst-kinetics-structural-grounding-2026-allow.txt`

---

## 1. Human CFTR Cryo-EM Structural Coordinates (`PDB 6MSM` vs. `PDB 5UAK`)

High-resolution cryo-electron microscopy (`cryo-EM`) reconstructions of full-length human `CFTR` (1,480 amino acids) anchor the conformational transition from the dephosphorylated, ATP-free inward-facing state (`rcsb-pdb-5uak`) to the phosphorylated, ATP-bound NBD-dimerized state (`rcsb-pdb-6msm`):

| Structural State & Catalog ID | PDB ID | Resolution ($\text{\AA}$) | Placed Residues (Chain A) | Unplaced Disordered Residues | Regulatory (`R`) Domain Coordinate Status |
| --- | --- | --- | --- | --- | --- |
| Dephosphorylated, ATP-free human `CFTR` (`rcsb-pdb-5uak`) | `5UAK` | 3.87 | 1,139 of 1,480 | 341 (`1480 - 1139 = 341`) | Atomic coordinates for the intrinsically disordered central portion of the `R` domain are `"not stated"` in the coordinate record |
| Phosphorylated, ATP-bound human `CFTR` `E1371Q` (`rcsb-pdb-6msm`) | `6MSM` | 3.2 | 1,181 of 1,480 | 299 (`1480 - 1181 = 299`) | Atomic coordinates for the intrinsically disordered central portion of the `R` domain are `"not stated"` in the coordinate record |
| Conformational Ordering Delta (`6MSM` minus `5UAK`) | — | — | +42 ordered residues (`1181 - 1139 = 42`) | -42 unplaced residues (`341 - 299 = 42`) | Reflects disengagement of the `R` domain from the intracellular vestibule and tightening of the `NBD1`-`NBD2` head-to-tail dimer |

Atomic coordinates for the intrinsically disordered central portion of the regulatory domain are not stated in the deposited PDB 6MSM or PDB 5UAK coordinate records because the flexible segments lack continuous cryo-EM map density.

---

## 2. Non-Equilibrium Irreversible Gating Cycle & Burst-Duration Coupling

Unlike classic ligand-gated channels that obey microscopic reversibility at thermal equilibrium, `CFTR` is an `ABC` transporter evolved into an ATP-driven chloride channel whose gating cycle is coupled to non-equilibrium nucleotide hydrolysis:

1. **Irreversible Gating Cycle & `E1371Q` Hydrolytic Arrest (`csanady-2006-jgp-irreversible-cycle`, DOI: `10.1085/jgp.200609558`):** Csanady L et al. 2006 demonstrated in the *Journal of General Physiology* that wild-type `CFTR` channel gating follows a non-equilibrium irreversible cycle driven by ATP binding at `NBD1` and `NBD2` and ATP hydrolysis at canonical composite site 2. In the hydrolysis-deficient `E1371Q` catalytic-glutamate mutant, blocking ATP hydrolysis at composite site 2 slows channel closure from open bursts by more than two orders of magnitude, prolonging mean open burst duration from approximately 350 milliseconds in wild-type `CFTR` to over 35,000 milliseconds (`35000 / 350 = 100`-fold prolongation).
2. **Stoichiometric 1:1 Burst-Hydrolysis Coupling (`csanady-2010-pnas-burst-coupling`, DOI: `10.1073/pnas.0911061107`):** Csanady L et al. 2010 demonstrated in *PNAS* strict stoichiometric 1:1 coupling between the `CFTR` open-burst state and the hydrolysis of a single ATP molecule at canonical composite site 2, showing that open burst duration distributions exhibit a peaked non-monotonic probability density with a deficit of very short openings.
3. **Single-Channel Conductance & Open Probability Benchmark (`csanady-2019-physrev-cftr-review`, DOI: `10.1152/physrev.00007.2018`):** Csanady L, Vergani P, and Gadsby DC 2019 synthesized in *Physiological Reviews* that wild-type human `CFTR` exhibits a single-channel chloride conductance of 8 to 10 pS and a maximal open probability of approximately 0.40 in excised inside-out membrane patches at saturating 2 mM MgATP after PKA phosphorylation.

---

## 3. Structural-Kinetic Synthesis

Combining the 100-fold open-burst prolongation in `E1371Q` (from 350 to 35,000 milliseconds) with the 1,181 ordered residues in `PDB 6MSM`, we infer that the `E1371Q` catalytic-glutamate substitution traps the `NBD1`-`NBD2` dimer in the pre-hydrolytic open-burst conformation required for high-resolution cryo-EM reconstruction. At the same time, because atomic coordinates for the intrinsically disordered central portion of the regulatory (`R`) domain are `"not stated"` in the deposited `PDB 6MSM` or `PDB 5UAK` coordinate records because the flexible segments lack continuous cryo-EM map density, any structural model of `R`-domain phosphorylation must explicitly distinguish experimentally placed residues (1,181 in `6MSM`; 1,139 in `5UAK`) from unplaced conformational ensembles.
