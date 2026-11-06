# CF Research Discovery Workbench (`tools/discovery/`)

> **Surface Organization (One-Primary-Axis Rule):** Inside `cf-research/tools/discovery/`, the primary navigation axis is **Scientific Discipline** (`structural_biology`, `rna_gene_editing`, `airway_delivery`, `microbiology_host`, `clinical_trials`, `ai_verification`), while the target personas are **CF Researchers** (`researchers`) and **Builders** (`builders`). Synthetic ODE/parameter-sweep modules are kept strictly in this researcher workbench and out of the `Patients & Families` hub on `senseiewok.ai`.

## Triage & Consolidation Summary

- **Total Verified Modules:** `103` (plus 4 infrastructure scripts and 165 unit test suites)
- **`keep` (Canonical Core Workbenches):** `85`
- **`merge` (Consolidate into Canonical Module):** `14`
- **`archive` (Exploratory / Superseded Parameter Sweeps):** `4`
- **Formal Proof & Stochastic Thermodynamics Bridge Suites (`formal_proof_thermodynamic_landscape.py`):** `59`

Verify this index, all module registrations, and the pure-stdlib unit test suites deterministically with:

```bash
python tools/discovery/discovery_council.py --audit-registry
python tools/discovery/discovery_council.py --check-index
python tools/discovery/discovery_council.py --falsify-all
python tools/discovery/discovery_council.py --run-stdlib-tests
```

## Structural Biology & CFTR Biophysics (`structural_biology` — 14 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `allosteric_potentiator_model.py` | `T-0221` | `T2` | `merge` | Merge into `gating_potentiator_energy_coupling.py` |
| `alpha_evolve_compound_het_ptc.py` | `T-0432` | `T2` | `keep` | Canonical workbench |
| `cftr_interdomain_allosteric_energetics.py` | `T-0447` | `T2` | `keep` | Canonical workbench |
| `cftr_membrane_insertion_thermodynamics.py` | `T-0453` | `T2` | `keep` | Canonical workbench |
| `cftr_sequence_scanner.py` | `T-0264` | `T2` | `keep` | Canonical workbench |
| `cftr_truncation_corrector_rescuability.py` | `T-0465` | `T2` | `keep` | Canonical workbench |
| `cftr_variant_catalog_concordance.py` | `T-0420` | `T2` | `keep` | Canonical workbench |
| `cryoem_truncation_model.py` | `T-0151` | `T2` | `keep` | Canonical workbench |
| `gating_potentiator_energy_coupling.py` | `T-0270` | `T2` | `keep` | Canonical workbench |
| `nbd1_thermodynamics.py` | `T-0199` | `T2` | `keep` | Canonical workbench |
| `pdi_redox_remodeling.py` | `T-0281` | `T2` | `archive` | Candidate for archive |
| `ptc_cryoem_allostery.py` | `T-0217` | `T2` | `merge` | Merge into `cryoem_truncation_model.py` |
| `slc26a9_tmem16a_bypass.py` | `T-0256` | `T2` | `keep` | Canonical workbench |
| `tmem16a_enac_synergy_model.py` | `T-0430` | `T2` | `keep` | Canonical workbench |

## RNA, Nonsense Readthrough & Gene Editing (`rna_gene_editing` — 25 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `alpha_evolve_cassette_shield.py` | `T-0433` | `T2` | `keep` | Canonical workbench |
| `alpha_evolve_splice_aso.py` | `T-0434` | `T2` | `keep` | Canonical workbench |
| `aso_splice_kinetics.py` | `T-0254` | `T2` | `keep` | Canonical workbench |
| `aso_splice_model.py` | `T-0183` | `T2` | `merge` | Merge into `aso_splice_kinetics.py` |
| `base_editing_cftr.py` | `T-0219` | `T2` | `merge` | Merge into `class1_prime_base_editing.py` |
| `bridge_rna_minipromoter_coevolve.py` | `T-0427` | `T2` | `keep` | Canonical workbench |
| `cftr_cryptic_splice_aso_kinetics.py` | `T-0470` | `T2` | `keep` | Canonical workbench |
| `cftr_gquadruplex_stability_model.py` | `T-0316` | `T2` | `archive` | Candidate for archive |
| `cftr_minigene_model.py` | `T-0207` | `T2` | `keep` | Canonical workbench |
| `class1_prime_base_editing.py` | `T-0255` | `T2` | `keep` | Canonical workbench |
| `dual_aav_intein_model.py` | `T-0206` | `T2` | `keep` | Canonical workbench |
| `intronic_splice_foundation.py` | `T-0269` | `T2` | `keep` | Canonical workbench |
| `minidna_evasion_model.py` | `T-0233` | `T2` | `archive` | Candidate for archive |
| `mrna_hypoxia_translation_dynamics.py` | `T-0452` | `T2` | `keep` | Canonical workbench |
| `mrna_poly_a_deadenylation_kinetics.py` | `T-0449` | `T2` | `keep` | Canonical workbench |
| `mrna_translation_velocity.py` | `T-0272` | `T2` | `keep` | Canonical workbench |
| `nmd_ptc_kinetics.py` | `T-0241` | `T2` | `merge` | Merge into `ptc_erf_depletion_readthrough.py` |
| `prime_editing_cftr.py` | `T-0184` | `T2` | `merge` | Merge into `class1_prime_base_editing.py` |
| `prime_editing_ptc_nmd_kinetics.py` | `T-0440` | `T2` | `keep` | Canonical workbench |
| `ptc_erf_depletion_readthrough.py` | `T-0401` | `T2` | `keep` | Canonical workbench |
| `split_intein_cftr_reconstitution.py` | `T-0393` | `T2` | `keep` | Canonical workbench |
| `suppressor_trna_model.py` | `T-0216` | `T2` | `keep` | Canonical workbench |
| `suppressor_trna_orthogonality.py` | `T-0394` | `T2` | `keep` | Canonical workbench |
| `u1_snrnp_decoy_kinetics.py` | `T-0292` | `T2` | `merge` | Merge into `u1_snrnp_decoy_thermodynamics.py` |
| `u1_snrnp_decoy_thermodynamics.py` | `T-0352` | `T2` | `keep` | Canonical workbench |

## Airway Delivery, Mucus Rheology & Alternative Ion Transport (`airway_delivery` — 51 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `aav_tropism_model.py` | `T-0208` | `T2` | `keep` | Canonical workbench |
| `aerosol_droplet_deposition_mechanics.py` | `T-0408` | `T2` | `keep` | Canonical workbench |
| `aerosol_hygroscopic_deposition_kinetics.py` | `T-0473` | `T2` | `keep` | Canonical workbench |
| `aerosol_mmad_mucus_penetration.py` | `T-0547` | `T2` | `keep` | Canonical workbench |
| `airway_basal_tight_junction_kinetics.py` | `T-0429` | `T2` | `keep` | Canonical workbench |
| `ali_airway_lnp_transfection_kinetics.py` | `T-0446` | `T2` | `keep` | Canonical workbench |
| `antimicrobial_peptide_acidification.py` | `T-0461` | `T2` | `keep` | Canonical workbench |
| `asl_osmotic_ciliary_transport.py` | `T-0459` | `T2` | `keep` | Canonical workbench |
| `asl_ph_antimicrobial_potency.py` | `T-0469` | `T2` | `keep` | Canonical workbench |
| `basal_stem_persistence.py` | `T-0220` | `T2` | `keep` | Canonical workbench |
| `basal_stem_quiescence_kinetics.py` | `T-0414` | `T2` | `keep` | Canonical workbench |
| `basal_stem_tropism_persistence.py` | `T-0498` | `T2` | `keep` | Canonical workbench |
| `biophysical_models.py` | `T-0154` | `T2` | `keep` | Canonical workbench |
| `ciliary_wave_clearance_model.py` | `T-0315` | `T2` | `keep` | Canonical workbench |
| `enac_knockdown_model.py` | `T-0200` | `T2` | `merge` | Merge into `enac_tmem16a_equilibrium.py` |
| `enac_proteolysis_asl_regulation.py` | `T-0417` | `T2` | `keep` | Canonical workbench |
| `enac_proteolysis_model.py` | `T-0245` | `T2` | `merge` | Merge into `enac_tmem16a_equilibrium.py` |
| `enac_tmem16a_equilibrium.py` | `T-0252` | `T2` | `keep` | Canonical workbench |
| `endocytic_pathway_partitioning.py` | `T-0561` | `T2` | `keep` | Canonical workbench |
| `endocytosis_mucociliary_competition.py` | `T-0416` | `T2` | `keep` | Canonical workbench |
| `endosomal_escape_kinetics.py` | `T-0557` | `T2` | `keep` | Canonical workbench |
| `epithelial_plasticity.py` | `T-0178` | `T2` | `keep` | Canonical workbench |
| `ev_basal_tropism.py` | `T-0234` | `T2` | `archive` | Candidate for archive |
| `glycocalyx_charge_barrier_model.py` | `T-0351` | `T2` | `merge` | Merge into `glycocalyx_permeabilization_model.py` |
| `glycocalyx_permeabilization_model.py` | `T-0291` | `T2` | `keep` | Canonical workbench |
| `lentiviral_basal_integration_persistence.py` | `T-0407` | `T2` | `keep` | Canonical workbench |
| `lnp_escape_dynamics.py` | `T-0285` | `T2` | `keep` | Canonical workbench |
| `lnp_mucus_shear_shedding_kinetics.py` | `T-0439` | `T2` | `keep` | Canonical workbench |
| `lnp_nuclear_entry_basal_persistence.py` | `T-0478` | `T2` | `keep` | Canonical workbench |
| `lnp_peg_shedding_mucus_kinetics.py` | `T-0463` | `T2` | `keep` | Canonical workbench |
| `lnp_stem_endosomal_escape_kinetics.py` | `T-0455` | `T2` | `keep` | Canonical workbench |
| `lnp_surfactant_corona_remodeling.py` | `T-0456` | `T2` | `keep` | Canonical workbench |
| `macrophage_evasion.py` | `T-0544` | `T2` | `keep` | Canonical workbench |
| `mechano_acoustic_delivery.py` | `T-0228` | `T2` | `keep` | Canonical workbench |
| `mpp_obstruction_scaling.py` | `T-0244` | `T2` | `keep` | Canonical workbench |
| `mrna_lnp_escape_model.py` | `T-0239` | `T2` | `merge` | Merge into `lnp_escape_dynamics.py` |
| `mrna_lnp_nebulizer_shear_stability.py` | `T-0481` | `T2` | `keep` | Canonical workbench |
| `mucociliary_clearance_viscoelasticity.py` | `T-0558` | `T2` | `keep` | Canonical workbench |
| `mucus_cleavage_kinetics.py` | `T-0253` | `T2` | `merge` | Merge into `mucus_rheology.py` |
| `mucus_penetrating_aso_kinetics.py` | `T-0497` | `T2` | `keep` | Canonical workbench |
| `mucus_rheology.py` | `T-0185` | `T2` | `keep` | Canonical workbench |
| `nebulizer_shear_stress_integrity.py` | `T-0413` | `T2` | `keep` | Canonical workbench |
| `periciliary_brush_compression.py` | `T-0562` | `T2` | `keep` | Canonical workbench |
| `surfactant_shielded_rnp_mucus_transport.py` | `T-0482` | `T2` | `keep` | Canonical workbench |
| `tmem16a_calcium_conductance_model.py` | `T-0392` | `T2` | `keep` | Canonical workbench |
| `tmem16a_conductance_model.py` | `T-0240` | `T2` | `merge` | Merge into `tmem16a_calcium_conductance_model.py` |
| `tmem16a_enac_asl_cbf_dynamics.py` | `T-0444` | `T2` | `keep` | Canonical workbench |
| `tmem16a_slc26a9_conductance_coupling.py` | `T-0476` | `T2` | `keep` | Canonical workbench |
| `transcytosis_model.py` | `T-0177` | `T2` | `keep` | Canonical workbench |
| `vector_redosing_kinetics.py` | `T-0198` | `T2` | `keep` | Canonical workbench |
| `vx522_inhaled_lnp_mucus_diffusion.py` | `T-0400` | `T2` | `keep` | Canonical workbench |

## Microbiology, Phage Therapy & Host Defense (`microbiology_host` — 7 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `alginate_lyase_phage_synergy.py` | `T-0472` | `T2` | `keep` | Canonical workbench |
| `macrophage_efferocytosis_model.py` | `T-0251` | `T2` | `merge` | Merge into `macrophage_evasion_kinetics.py` |
| `macrophage_evasion_kinetics.py` | `T-0271` | `T2` | `keep` | Canonical workbench |
| `phage_biofilm_kinetics.py` | `T-0197` | `T2` | `keep` | Canonical workbench |
| `phage_biofilm_penetration_coevolution.py` | `T-0445` | `T2` | `keep` | Canonical workbench |
| `phage_depolymerase_biofilm_synergy.py` | `T-0402` | `T2` | `keep` | Canonical workbench |
| `phage_omko1_tradeoff.py` | `T-0246` | `T2` | `keep` | Canonical workbench |

## Clinical Trials, Registries & Theratyping (`clinical_trials` — 3 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `counterfactual_trial_forensics.py` | `T-0428` | `T2` | `keep` | Canonical workbench |
| `modifier_interactome.py` | `T-0186` | `T2` | `keep` | Canonical workbench |
| `rare_allele_matrix.py` | `T-0201` | `T2` | `keep` | Canonical workbench |

## AI Verification & Scientific Benchmarking (`ai_verification` — 3 modules)

| Module | Task | Tier | Triage | Consolidation Target / Notes |
| --- | --- | --- | --- | --- |
| `alpha_evolve_pareto_regimen.py` | `T-0426` | `T2` | `keep` | Canonical workbench |
| `formal_proof_thermodynamic_landscape.py` | `T-0479` | `T2` | `keep` | Canonical workbench |
| `science_model_benchmark.py` | `T-0265` | `T2` | `keep` | Canonical workbench |

## Infrastructure, Cross-Domain Mining & Symbolic Verification (`4` scripts + `2` cross-cutting suites)

| Script / Suite | Role | Companion Verification |
| --- | --- | --- |
| `discovery_council.py` | 3-phase Curative Discovery Council runner, dynamic module registry auditor, and AH-10 falsification harness | `test_discovery_council.py` |
| `cross_domain_miner.py` | Cross-domain delivery vector feasibility (BFI), transcytosis, and basal stem-cell trajectory miner | `test_cross_domain_miner.py` |
| `export_delivery_web_data.py` | Deterministic TSV/JSON snapshot exporter for `data/delivery_discovery_snapshot.json` | `test_export_delivery_web_data.py` |
| `epistemic_rigor.py` | Shared epistemic contract, sensitivity perturbation, and negative-control verification helpers | `test_epistemic_rigor_and_falsification.py` |
| `test_symbolic_invariants.py` | Pure-stdlib complex-step autodiff (`Im(f(x+ih))/h`), exact integrals, and Buckingham Pi dimensional homogeneity (`T-0505`) | Self-contained (`unittest`) |
| `data/canonical_math_bridge_snapshot.json` | Locked Lean 4 / formal proof structural snapshot (`fermats_last_theorem`, `formal_math_zeta23`, `formal_math_3sum_apsp`) | `test_formal_proof_thermodynamic_landscape.py` |

## Formal Proof & Stochastic Thermodynamics Bridge Suites (`formal_proof_thermodynamic_landscape.py` — 59 suites)

| Test Suite | Round | Task | Verified Theorem & CF Biophysical Bridge |
| --- | --- | --- | --- |
| `test_barato_hartich_seifert_dissipation.py` | `Round 50` | `T-0605` | Unit tests for Round 50 (T-0605): Barato-Hartich-Seifert Information Thermodynamics & Dissipation Bounds in CFTR Gating in cf-research |
| `test_barato_seifert_brownian_clock_precision.py` | `Round 82` | `T-0675` | Unit tests for Round 82 (T-0675): Barato-Seifert Cost of Precision & Biomolecular Clock Efficiency in CFTR Catalytic Assemblies in cf-research |
| `test_barato_seifert_universal_tur.py` | `Round 74` | `T-0659` | Unit tests for Round 74 (T-0659): Barato-Seifert Universal Thermodynamic Uncertainty Relation for Molecular Motors in CFTR Catalytic Assemblies in cf-research |
| `test_brandner_hanazato_saito_ballistic_tur.py` | `Round 85` | `T-0682` | Unit tests for Round 85 (T-0682): Brandner-Hanazato-Saito Universal Thermodynamic Bounds on Ballistic Transport & Channel Precision in cf-research |
| `test_current_fluctuation_tur.py` | `Round 45` | `T-0595` | Unit tests for Round 45 (T-0595): Seifert-Gingrich Current Fluctuation Theorems & Generalized TUR for CFTR Channel Kinetics in cf-research |
| `test_dechant_sasa_bounds.py` | `Round 44` | `T-0593` | Unit tests for Round 44 (T-0593): Dechant-Sasa Wasserstein Fluctuation Dissipation Bounds for CFTR Gating in cf-research |
| `test_dechant_sasa_fluctuation_response.py` | `Round 58` | `T-0620` | Unit tests for Round 58 (T-0620): Dechant-Sasa Fluctuation-Response Inequality and Out-of-Equilibrium Susceptibility in CFTR Gating in cf-research |
| `test_dechant_sasa_multidimensional_tur.py` | `Round 75` | `T-0661` | Unit tests for Round 75 (T-0661): Dechant-Sasa Multidimensional Thermodynamic Uncertainty Relations & Fluctuation-Response Bounds in CFTR Catalytic Assemblies in cf-research |
| `test_dechant_sasa_speed_limits.py` | `Round 52` | `T-0607` | Unit tests for Round 52 (T-0607): Dechant-Sasa Thermodynamic Speed Limits & Geometric Bounds in CFTR Gating in cf-research |
| `test_dechant_sasa_wasserstein_bounds.py` | `Round 65` | `T-0635` | Unit tests for Round 65 (T-0635): Dechant-Sasa Wasserstein Distance & Geometric Fluctuation Bounds in CFTR Macromolecular Transition Manifolds in cf-research |
| `test_esposito_parrondo_sagawa_information_efficiency.py` | `Round 67` | `T-0639` | Unit tests for Round 67 (T-0639): Esposito-Parrondo-Sagawa Second-Law Information Efficiency & Allosteric Feedback Bounds in CFTR Macromolecular Gating Cascades in cf-research |
| `test_esposito_van_den_broeck_entropy.py` | `Round 54` | `T-0612` | Unit tests for Round 54 (T-0612): Esposito-Van den Broeck Three-Way Entropy Production Decomposition & Housekeeping Dissipation in CFTR Gating in cf-research |
| `test_falasco_esposito_dissipation_time.py` | `Round 57` | `T-0619` | Unit tests for Round 57 (T-0619): Falasco-Esposito Dissipation-Time Uncertainty Relation in CFTR Non-Equilibrium Transitions in cf-research |
| `test_falasco_esposito_dissipation_time_uncertainty.py` | `Round 70` | `T-0646` | Unit tests for Round 70 (T-0646): Falasco-Esposito Dissipation-Time Uncertainty Relations & Gating Precision Trade-offs in Driven CFTR Catalytic Assemblies in cf-research |
| `test_freidlin_wentzell_quasipotential.py` | `Round 30` | `T-0574` | Unit tests for Round 30 (T-0574): Freidlin-Wentzell Large Deviation Theory & Quasipotential Action for Rare CFTR Opening Events |
| `test_garrahan_frenetic_first_passage_bound.py` | `Round 83` | `T-0677` | Unit tests for Round 83 (T-0677): Garrahan Dynamical Activity First-Passage Bounds & Frenetic Gating Regularity in CFTR Catalytic Assemblies in cf-research |
| `test_geometric_gauge_holonomy.py` | `Round 37` | `T-0581` | Unit tests for Round 37 (T-0581): Geometrical Thermodynamic Gauge Theory & Berry Holonomy for CFTR Catalytic Cycles in cf-research |
| `test_gingrich_horowitz_first_passage_tur.py` | `Round 81` | `T-0673` | Unit tests for Round 81 (T-0673): Gingrich-Horowitz First Arrival Time Bounds & Dwell-Time TUR in CFTR Catalytic Assemblies in cf-research |
| `test_gingrich_horowitz_thermodynamic_uncertainty.py` | `Round 72` | `T-0655` | Unit tests for Round 72 (T-0655): Gingrich-Horowitz Finite-Time Thermodynamic Uncertainty Relations & Current Precision Bounds in CFTR Catalytic Assemblies in cf-research |
| `test_gingrich_rotskoff_dissipation_equality.py` | `Round 79` | `T-0669` | Unit tests for Round 79 (T-0669): Gingrich-Rotskoff-Vaikuntanathan-Horowitz Dissipation-Uncertainty Equality in CFTR Catalytic Assemblies in cf-research |
| `test_glansdorff_prigogine_evolution.py` | `Round 26` | `T-0553` | Unit tests for Round 26 (T-0553): Glansdorff-Prigogine Universal Evolution Criterion & Non-Linear Gating Dynamics in CFTR |
| `test_graham_haken_potential.py` | `Round 31` | `T-0575` | Unit tests for Round 31 (T-0575): Graham-Haken Non-Equilibrium Generalized Potential for Multi-Stable CFTR Gating |
| `test_hartich_barato_seifert_bipartite.py` | `Round 47` | `T-0599` | Unit tests for Round 47 (T-0599): Hartich-Barato-Seifert Bipartite Information Thermodynamics & Subsystem Entropy Bounds in cf-research |
| `test_hatano_sasa_steady_state_thermodynamics.py` | `Round 55` | `T-0615` | Unit tests for Round 55 (T-0615): Hatano-Sasa Steady-State Thermodynamics & Dissipated Work Bounds in CFTR Gating Transitions in cf-research |
| `test_horowitz_parrondo_feedback_control.py` | `Round 51` | `T-0606` | Unit tests for Round 51 (T-0606): Horowitz-Parrondo Information Flow & Feedback Control in CFTR Catalytic Assemblies in cf-research |
| `test_knst_nonequilibrium_work.py` | `Round 35` | `T-0579` | Unit tests for Round 35 (T-0579): Komatsu-Nakagawa-Sasa-Tasaki (KNST) Non-Equilibrium Work & Free Energy for CFTR Gating Transitions |
| `test_koyuk_seifert_biochemical_clock_speed.py` | `Round 69` | `T-0643` | Unit tests for Round 69 (T-0643): Koyuk-Seifert Thermodynamic Speed Bounds on Biochemical Clocks & Cyclic CFTR Macromolecular Oscillations in cf-research |
| `test_koyuk_seifert_cyclic_speed_dissipation.py` | `Round 56` | `T-0616` | Unit tests for Round 56 (T-0616): Koyuk-Seifert General Bound on Speed, Precision, and Dissipation in CFTR Biochemical Gating Cycles in cf-research |
| `test_koyuk_seifert_operationally_accessible_bounds.py` | `Round 76` | `T-0663` | Unit tests for Round 76 (T-0663): Koyuk-Seifert Operationally Accessible Bounds on Fluctuations and Dissipation in CFTR Catalytic Assemblies in cf-research |
| `test_koyuk_seifert_periodic_tur.py` | `Round 84` | `T-0681` | Unit tests for Round 84 (T-0681): Koyuk-Seifert Periodically Driven Thermodynamic Uncertainty Relations in Cyclic CFTR Catalytic Assemblies in cf-research |
| `test_kpv_relative_entropy.py` | `Round 43` | `T-0591` | Unit tests for Round 43 (T-0591): Kawai-Parrondo-Van den Broeck (KPV) Dissipation-Information Relative Entropy Equality for CFTR Trajectories in cf-research |
| `test_lebowitz_spohn_action.py` | `Round 59` | `T-0623` | Unit tests for Round 59 (T-0623): Lebowitz-Spohn Stochastic Action Functional and Gallavotti-Cohen Symmetry in CFTR Gating Transitions in cf-research |
| `test_li_horowitz_dissipation_quantification.py` | `Round 80` | `T-0671` | Unit tests for Round 80 (T-0671): Li-Horowitz-Gingrich-Fakhri Dissipation Quantification from Fluctuations in CFTR Catalytic Assemblies in cf-research |
| `test_maes_netocny_response.py` | `Round 33` | `T-0577` | Unit tests for Round 33 (T-0577): Maes-Netočný Non-Equilibrium Response Formula & Excess Entropy Production for Rescued CFTR Kinetics |
| `test_manikandan_gupta_krishnamurthy_bounds.py` | `Round 63` | `T-0631` | Unit tests for Round 63 (T-0631): Manikandan-Gupta-Krishnamurthy Operationally Accessible Bounds on Dissipation & Current Fluctuations in CFTR Channel Networks in cf-research |
| `test_onsager_electroosmotic_transport.py` | `Round 24` | `T-0546` | Unit tests for Round 24 (T-0546): Onsager Reciprocal Relations & Cross-Coupled Electro-Osmotic Transport Across Airway Surface Liquid |
| `test_onsager_machlup_action.py` | `Round 29` | `T-0567` | Unit tests for Round 29 (T-0567): Stochastic Path Integral & Onsager-Machlup Action for Active CFTR Dimerization |
| `test_parrondo_horowitz_sagawa_information.py` | `Round 48` | `T-0601` | Unit tests for Round 48 (T-0601): Parrondo-Horowitz-Sagawa Information Thermodynamics & Maxwell Demon Feedback Cooling in CFTR Channel Assemblies in cf-research |
| `test_pietzonka_barato_seifert_current_bounds.py` | `Round 78` | `T-0667` | Unit tests for Round 78 (T-0667): Pietzonka-Barato-Seifert Universal Finite-Time Current Bounds in CFTR Catalytic Assemblies in cf-research |
| `test_pietzonka_seifert_bounds.py` | `Round 46` | `T-0597` | Unit tests for Round 46 (T-0597): Pietzonka-Seifert Finite-Time Dissipation Bounds & Molecular Machine Efficiency in cf-research |
| `test_pietzonka_seifert_power_efficiency_constancy.py` | `Round 73` | `T-0657` | Unit tests for Round 73 (T-0657): Pietzonka-Seifert Universal Power-Efficiency-Constancy Trade-Off & Molecular Engine Bounds in CFTR Catalytic Assemblies in cf-research |
| `test_polettini_esposito_cycle_decomposition.py` | `Round 61` | `T-0627` | Unit tests for Round 61 (T-0627): Polettini-Esposito Non-Convex Relative Entropy and Cycle Dissipation Decomposition in CFTR Gating Cycles in cf-research |
| `test_polettini_geometric_gauge_decomposition.py` | `Round 66` | `T-0636` | Unit tests for Round 66 (T-0636): Polettini Geometric Gauge Decomposition & Non-Convex Relative Entropy in CFTR Macromolecular Transition Networks in cf-research |
| `test_prigogine_minimum_entropy_production.py` | `Round 25` | `T-0551` | Unit tests for Round 25 (T-0551): Prigogine Minimum Entropy Production Principle & Stationary NESS Stability in CFTR Gating |
| `test_proesmans_discrete_time_tur.py` | `Round 86` | `T-0685` | Unit tests for Round 86 (T-0685): Proesmans-Van den Broeck Discrete-Time TUR & Stroboscopic CFTR Gating Bounds in cf-research |
| `test_prost_joanny_parrondo_fdt.py` | `Round 41` | `T-0587` | Unit tests for Round 41 (T-0587): Prost-Joanny-Parrondo Generalized Fluctuation-Dissipation Theorem for Molecular Motors & CFTR Active Gating in cf-research |
| `test_rao_esposito_conservation_laws.py` | `Round 62` | `T-0629` | Unit tests for Round 62 (T-0629): Rao-Esposito Conservation Laws and Topological Decomposition of Dissipation in CFTR Channel Networks in cf-research |
| `test_sagawa_ueda_fluctuation_efficiency.py` | `Round 49` | `T-0603` | Unit tests for Round 49 (T-0603): Sagawa-Ueda Generalized Fluctuation Theorem & Information Efficacy in Allosteric Drug Stabilization in cf-research |
| `test_sasa_tasaki_equality.py` | `Round 39` | `T-0583` | Unit tests for Round 39 (T-0583): Non-Equilibrium Steady-State Transition Bounds & Sasa-Tasaki Equality in CFTR Catalytic Assemblies in cf-research |
| `test_schnakenberg_cycle_kinetics.py` | `Round 28` | `T-0565` | Unit tests for Round 28 (T-0565): Schnakenberg Cycle Kinetics, Non-Equilibrium Cycle Flux Splitting & Detailed Balance Breaking in CFTR Gating |
| `test_seifert_trajectory_entropy_production.py` | `Round 40` | `T-0585` | Unit tests for Round 40 (T-0585): Seifert Trajectory Stochastic Entropy Production & Fluctuation Theorem in CFTR Catalytic Cycles in cf-research |
| `test_shiraishi_funo_saito_speed_limit.py` | `Round 68` | `T-0641` | Unit tests for Round 68 (T-0641): Shiraishi-Funo-Saito Thermodynamic Speed Limits & Operational Gating Bounds in CFTR Macromolecular Transition Dynamics in cf-research |
| `test_shiraishi_funo_saito_speed_limits.py` | `Round 53` | `T-0611` | Unit tests for Round 53 (T-0611): Shiraishi-Funo-Saito Universal Thermodynamic Speed Limits & Operational Cost Bounds in CFTR Channel Gating in cf-research |
| `test_shiraishi_saito_tasaki_power_efficiency.py` | `Round 87` | `T-0686` | Unit tests for Round 87 (T-0686): Shiraishi-Saito-Tasaki Power-Efficiency Trade-Off & CFTR NBD1-NBD2 Engine Bounds in cf-research |
| `test_speck_seifert_fdt.py` | `Round 42` | `T-0589` | Unit tests for Round 42 (T-0589): Speck-Seifert Non-Equilibrium Fluctuation-Dissipation Theorem & Velocity Field Coupling for CFTR Active Gating in cf-research |
| `test_speck_seifert_housekeeping.py` | `Round 60` | `T-0624` | Unit tests for Round 60 (T-0624): Speck-Seifert Housekeeping Heat Fluctuation Theorem and Non-Equilibrium Dissipation in CFTR Gating Cycles in cf-research |
| `test_van_vu_hasegawa_geometrical_irreversibility.py` | `Round 77` | `T-0665` | Unit tests for Round 77 (T-0665): Van Vu-Hasegawa Geometrical Bounds of Irreversibility and Dynamical Activity in CFTR Catalytic Assemblies in cf-research |
| `test_van_vu_hasegawa_thermodynamic_uncertainty.py` | `Round 71` | `T-0650` | Unit tests for Round 71 (T-0650): Van Vu - Hasegawa Thermodynamic Uncertainty Relations & Dynamical Activity Bounds in Arbitrarily Driven CFTR Catalytic Assemblies in cf-research |
| `test_wachtel_rao_esposito_coarse_graining.py` | `Round 64` | `T-0633` | Unit tests for Round 64 (T-0633): Wachtel-Rao-Esposito Thermodynamically Consistent Coarse Graining & Cycle Affinity Invariants in CFTR Conformational Transition Networks in cf-research |
