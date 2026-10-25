# CF Research Discovery Workbench (`tools/discovery/`)

> **Surface Organization (One-Primary-Axis Rule):** Inside `cf-research/tools/discovery/`, the primary navigation axis is **Scientific Discipline** (`structural_biology`, `rna_gene_editing`, `airway_delivery`, `microbiology_host`, `clinical_trials`, `ai_verification`), while the target personas are **CF Researchers** (`researchers`) and **Builders** (`builders`). Synthetic ODE/parameter-sweep modules are kept strictly in this researcher workbench and out of the `Patients & Families` hub on `senseiewok.ai`.

## Triage & Consolidation Summary

- **Total Verified Modules:** `57` (plus 3 infrastructure scripts and 58 unit test suites)
- **`keep` (Canonical Core Workbenches):** `39`
- **`merge` (Consolidate into Canonical Module):** `14`
- **`archive` (Exploratory / Superseded Parameter Sweeps):** `4`

Verify this index and all module registrations deterministically with:

```bash
python tools/discovery/discovery_council.py --audit-registry
python tools/discovery/discovery_council.py --check-index
```

## Structural Biology & CFTR Biophysics (`structural_biology` — 7 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `allosteric_potentiator_model.py` | `T-0221` | `T2` | `merge` | Merge into `gating_potentiator_energy_coupling.py` |
| `cftr_sequence_scanner.py` | `T-0264` | `T2` | `keep` | Canonical workbench |
| `cryoem_truncation_model.py` | `T-0151` | `T2` | `keep` | Canonical workbench |
| `gating_potentiator_energy_coupling.py` | `T-0270` | `T2` | `keep` | Canonical workbench |
| `nbd1_thermodynamics.py` | `T-0199` | `T2` | `keep` | Canonical workbench |
| `pdi_redox_remodeling.py` | `T-0281` | `T2` | `archive` | Candidate for archive |
| `ptc_cryoem_allostery.py` | `T-0217` | `T2` | `merge` | Merge into `cryoem_truncation_model.py` |

## RNA, Nonsense Readthrough & Gene Editing (`rna_gene_editing` — 17 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `aso_splice_kinetics.py` | `T-0254` | `T2` | `keep` | Canonical workbench |
| `aso_splice_model.py` | `T-0183` | `T2` | `merge` | Merge into `aso_splice_kinetics.py` |
| `base_editing_cftr.py` | `T-0219` | `T2` | `merge` | Merge into `class1_prime_base_editing.py` |
| `cftr_gquadruplex_stability_model.py` | `T-0316` | `T2` | `archive` | Candidate for archive |
| `cftr_minigene_model.py` | `T-0207` | `T2` | `keep` | Canonical workbench |
| `class1_prime_base_editing.py` | `T-0255` | `T2` | `keep` | Canonical workbench |
| `dual_aav_intein_model.py` | `T-0206` | `T2` | `keep` | Canonical workbench |
| `intronic_splice_foundation.py` | `T-0269` | `T2` | `keep` | Canonical workbench |
| `minidna_evasion_model.py` | `T-0233` | `T2` | `archive` | Candidate for archive |
| `mrna_translation_velocity.py` | `T-0272` | `T2` | `keep` | Canonical workbench |
| `nmd_ptc_kinetics.py` | `T-0241` | `T2` | `merge` | Merge into `ptc_erf_depletion_readthrough.py` |
| `prime_editing_cftr.py` | `T-0184` | `T2` | `merge` | Merge into `class1_prime_base_editing.py` |
| `ptc_erf_depletion_readthrough.py` | `T-0401` | `T2` | `keep` | Canonical workbench |
| `split_intein_cftr_reconstitution.py` | `T-0393` | `T2` | `keep` | Canonical workbench |
| `suppressor_trna_model.py` | `T-0216` | `T2` | `keep` | Canonical workbench |
| `u1_snrnp_decoy_kinetics.py` | `T-0292` | `T2` | `merge` | Merge into `u1_snrnp_decoy_thermodynamics.py` |
| `u1_snrnp_decoy_thermodynamics.py` | `T-0352` | `T2` | `keep` | Canonical workbench |

## Airway Delivery, Mucus Rheology & Alternative Ion Transport (`airway_delivery` — 25 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `aav_tropism_model.py` | `T-0208` | `T2` | `keep` | Canonical workbench |
| `aerosol_droplet_deposition_mechanics.py` | `T-0408` | `T2` | `keep` | Canonical workbench |
| `basal_stem_persistence.py` | `T-0220` | `T2` | `keep` | Canonical workbench |
| `basal_stem_quiescence_kinetics.py` | `T-0414` | `T2` | `keep` | Canonical workbench |
| `biophysical_models.py` | `T-0154` | `T2` | `keep` | Canonical workbench |
| `ciliary_wave_clearance_model.py` | `T-0315` | `T2` | `keep` | Canonical workbench |
| `enac_knockdown_model.py` | `T-0200` | `T2` | `merge` | Merge into `enac_tmem16a_equilibrium.py` |
| `enac_proteolysis_asl_regulation.py` | `T-0417` | `T2` | `keep` | Canonical workbench |
| `enac_proteolysis_model.py` | `T-0245` | `T2` | `merge` | Merge into `enac_tmem16a_equilibrium.py` |
| `enac_tmem16a_equilibrium.py` | `T-0252` | `T2` | `keep` | Canonical workbench |
| `endocytosis_mucociliary_competition.py` | `T-0416` | `T2` | `keep` | Canonical workbench |
| `epithelial_plasticity.py` | `T-0178` | `T2` | `keep` | Canonical workbench |
| `ev_basal_tropism.py` | `T-0234` | `T2` | `archive` | Candidate for archive |
| `glycocalyx_charge_barrier_model.py` | `T-0351` | `T2` | `merge` | Merge into `glycocalyx_permeabilization_model.py` |
| `glycocalyx_permeabilization_model.py` | `T-0291` | `T2` | `keep` | Canonical workbench |
| `lnp_escape_dynamics.py` | `T-0285` | `T2` | `keep` | Canonical workbench |
| `mechano_acoustic_delivery.py` | `T-0228` | `T2` | `keep` | Canonical workbench |
| `mrna_lnp_escape_model.py` | `T-0239` | `T2` | `merge` | Merge into `lnp_escape_dynamics.py` |
| `mucus_cleavage_kinetics.py` | `T-0253` | `T2` | `merge` | Merge into `mucus_rheology.py` |
| `mucus_rheology.py` | `T-0185` | `T2` | `keep` | Canonical workbench |
| `tmem16a_calcium_conductance_model.py` | `T-0392` | `T2` | `keep` | Canonical workbench |
| `tmem16a_conductance_model.py` | `T-0240` | `T2` | `merge` | Merge into `tmem16a_calcium_conductance_model.py` |
| `transcytosis_model.py` | `T-0177` | `T2` | `keep` | Canonical workbench |
| `vector_redosing_kinetics.py` | `T-0198` | `T2` | `keep` | Canonical workbench |
| `vx522_inhaled_lnp_mucus_diffusion.py` | `T-0400` | `T2` | `keep` | Canonical workbench |

## Microbiology, Phage Therapy & Host Defense (`microbiology_host` — 5 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `macrophage_efferocytosis_model.py` | `T-0251` | `T2` | `merge` | Merge into `macrophage_evasion_kinetics.py` |
| `macrophage_evasion_kinetics.py` | `T-0271` | `T2` | `keep` | Canonical workbench |
| `phage_biofilm_kinetics.py` | `T-0197` | `T2` | `keep` | Canonical workbench |
| `phage_depolymerase_biofilm_synergy.py` | `T-0402` | `T2` | `keep` | Canonical workbench |
| `phage_omko1_tradeoff.py` | `T-0246` | `T2` | `keep` | Canonical workbench |

## Clinical Trials, Registries & Theratyping (`clinical_trials` — 2 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `modifier_interactome.py` | `T-0186` | `T2` | `keep` | Canonical workbench |
| `rare_allele_matrix.py` | `T-0201` | `T2` | `keep` | Canonical workbench |

## AI Verification & Scientific Benchmarking (`ai_verification` — 1 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `science_model_benchmark.py` | `T-0265` | `T2` | `keep` | Canonical workbench |
