#!/usr/bin/env python3
"""Unit test suite for Airway Epithelial Glycocalyx Charge Barrier Model (T-0351).

Acceptance Criteria:
  pytest tools/discovery/test_glycocalyx_charge_barrier_model.py -v; python tools/discovery/glycocalyx_charge_barrier_model.py --check
"""
import math
import sys
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from glycocalyx_charge_barrier_model import (
    ASL_IONIC_STRENGTH_MM,
    BASELINE_PORE_RADIUS_NM,
    BASELINE_SURFACE_CHARGE_DENSITY,
    NLM_REGISTRY_DISCLAIMER,
    SIALIDASE_PROFILES,
    calculate_debye_length,
    calculate_surface_potential,
    evaluate_coupled_diffusivity,
    evaluate_electrostatic_double_layer,
    evaluate_steric_hindrance,
    export_complete_glycocalyx_manifest,
    run_checks,
    simulate_enzymatic_desialylation,
)


class TestGlycocalyxChargeBarrierModel(unittest.TestCase):
    """Deterministic biophysical tests for T-0351."""

    def test_debye_length_in_airway_surface_liquid(self):
        debye = calculate_debye_length(ASL_IONIC_STRENGTH_MM, 310.15)
        self.assertAlmostEqual(debye, 0.779, delta=0.05)

    def test_debye_length_ionic_strength_dependence(self):
        # Higher ionic strength -> shorter Debye screening length
        d_150 = calculate_debye_length(150.0, 310.15)
        d_300 = calculate_debye_length(300.0, 310.15)
        self.assertLess(d_300, d_150)
        self.assertAlmostEqual(d_300, d_150 / math.sqrt(2.0), delta=0.01)

    def test_surface_potential_grahame_calculation(self):
        psi_0 = calculate_surface_potential(BASELINE_SURFACE_CHARGE_DENSITY, 150.0, 310.15)
        # Baseline charge density (-0.035 C/m^2) should produce negative potential ~ -38 mV
        self.assertLess(psi_0, -20.0)
        self.assertGreater(psi_0, -50.0)

    def test_electrostatic_potential_decay(self):
        edl = evaluate_electrostatic_double_layer()
        self.assertLess(edl.surface_potential_mV, 0.0)
        # Potential magnitude decays exponentially with distance
        self.assertGreater(abs(edl.surface_potential_mV), abs(edl.potential_at_1nm_mV))
        self.assertGreater(abs(edl.potential_at_1nm_mV), abs(edl.potential_at_2nm_mV))
        self.assertGreater(abs(edl.potential_at_2nm_mV), abs(edl.potential_at_5nm_mV))

    def test_steric_partitioning_renkin_ogston(self):
        # r_p = 12.5 nm in 25 nm pore -> lambda = 0.5
        steric_12_5 = evaluate_steric_hindrance(12.5, pore_radius_nm=25.0)
        self.assertEqual(steric_12_5.pore_size_ratio_lambda, 0.5)
        self.assertAlmostEqual(steric_12_5.ogston_partition_phi, 0.25, delta=1e-3)
        self.assertGreater(steric_12_5.effective_steric_partition_k, 0.0)

        # Particle exceeding pore radius must be completely excluded
        steric_exceed = evaluate_steric_hindrance(30.0, pore_radius_nm=25.0)
        self.assertEqual(steric_exceed.effective_steric_partition_k, 0.0)

    def test_coupled_diffusivity_electrostatic_hindrance(self):
        # Neutral vs anionic particle of same size (r = 12.5 nm)
        diff_anionic = evaluate_coupled_diffusivity("Anionic_NP", 12.5, -20.0)
        diff_neutral = evaluate_coupled_diffusivity("Neutral_NP", 12.5, 0.0)

        # Electrostatic repulsion should reduce diffusivity of anionic particle
        self.assertLess(diff_anionic.k_elec, diff_neutral.k_elec)
        self.assertLess(diff_anionic.d_eff_glycocalyx_m2_s, diff_neutral.d_eff_glycocalyx_m2_s)

    def test_bacterial_neuraminidase_desialylation_kinetics(self):
        res = simulate_enzymatic_desialylation("bacterial_neuraminidase", enzyme_conc_uM=0.50, duration_minutes=30.0)
        self.assertGreaterEqual(res.desialylation_fraction, 0.80)
        self.assertGreaterEqual(res.diffusivity_fold_enhancement, 5.0)
        self.assertLess(res.t_80_percent_reduction_seconds, 600.0)  # less than 10 minutes

    def test_human_neuraminidase_isoforms(self):
        for enzyme in ["recombinant_human_neu1", "recombinant_human_neu2", "recombinant_human_neu3"]:
            res = simulate_enzymatic_desialylation(enzyme, enzyme_conc_uM=1.0, duration_minutes=60.0)
            self.assertGreater(res.desialylation_fraction, 0.50)
            self.assertGreater(res.diffusivity_fold_enhancement, 2.0)

    def test_manifest_grounding_and_provenance(self):
        manifest = export_complete_glycocalyx_manifest()
        self.assertEqual(manifest["model_id"], "T-0351")
        self.assertEqual(manifest["answer_tier"], "T2")
        self.assertEqual(manifest["track"], "Track A")
        self.assertIn("NCT05668741", manifest["clinical_trial_grounding"])
        self.assertIn("10.1126/science.1223012", manifest["primary_literature_dois"])
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", manifest["nlm_registry_disclaimer"])

    def test_self_check_runner(self):
        code = run_checks()
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
