#!/usr/bin/env python3
"""Deterministic verification test suite for modifier_interactome.py (T-0186).

Dual compatibility:
- pytest tools/discovery/test_modifier_interactome.py
- python -m unittest discover -s tools/discovery -p "test_*.py"
"""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import modifier_interactome as mi


class TestEpithelialIonFluxAndCurrents(unittest.TestCase):
    """Verifies apical epithelial ion flux balances and short-circuit current shifts."""

    def test_untreated_cf_absorptive_flux(self):
        res = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        self.assertEqual(res.currents.secretory_flux_direction, "ABSORPTIVE")
        self.assertEqual(res.currents.delta_isc_ua_cm2, 0.0)
        self.assertEqual(res.currents.net_short_circuit_current_isc_ua_cm2, -45.0)
        self.assertGreater(res.currents.i_na_ua_cm2, 40.0)

    def test_flux_direction_transition_to_secretory(self):
        res_dual = mi.simulate_modifier_regimen("dual_etd001_etx001", genotype_slc26a9="CT")
        self.assertEqual(res_dual.currents.secretory_flux_direction, "SECRETORY")
        self.assertGreaterEqual(res_dual.currents.delta_isc_ua_cm2, 30.0)

    def test_enac_blockade_reduces_na_current(self):
        res_base = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        res_etd = mi.simulate_modifier_regimen("monotherapy_etd001", genotype_slc26a9="CT")
        self.assertLess(res_etd.currents.i_na_ua_cm2, res_base.currents.i_na_ua_cm2 * 0.30)


class TestASLHeightRestoration(unittest.TestCase):
    """Verifies steady-state airway surface liquid height h_ASL dynamics and therapeutic rescue."""

    def test_untreated_cf_collapsed_asl(self):
        res = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        self.assertLess(res.hydration.steady_state_asl_height_um, 3.50)
        self.assertEqual(res.hydration.ciliary_beating_status, "CILIA_COMPRESSED_IMMOBILIZED")
        self.assertEqual(res.hydration.hydration_rescue_competence, "COLLAPSED")
        self.assertEqual(res.verdict, "INSUFFICIENT_ASL_COLLAPSED")

    def test_dual_combination_rehydration(self):
        # ETD001 + ETX001 combination restores ASL >= 7.0 um (physiological target)
        res_dual = mi.simulate_modifier_regimen("dual_etd001_etx001", genotype_slc26a9="CT")
        self.assertGreaterEqual(res_dual.hydration.steady_state_asl_height_um, mi.PHYSIOLOGICAL_ASL_TARGET_UM)
        self.assertEqual(res_dual.hydration.ciliary_beating_status, "UNIMPEDED_RAPID_CLEARANCE")
        self.assertEqual(res_dual.hydration.hydration_rescue_competence, "FULL_PHYSIOLOGICAL_RESCUE")
        self.assertEqual(res_dual.verdict, "CURATIVE_REHYDRATION_SUCCESS")

    def test_triple_bypass_superior_hydration(self):
        # Triple combination achieves even higher ASL height (> 8.0 um)
        res_triple = mi.simulate_modifier_regimen("triple_bypass", genotype_slc26a9="CT")
        self.assertGreater(res_triple.hydration.steady_state_asl_height_um, 8.0)
        self.assertEqual(res_triple.verdict, "CURATIVE_REHYDRATION_SUCCESS")

    def test_monotherapy_partial_rescue_bounds(self):
        # Monotherapies should improve ASL over untreated CF but stay below the 7.0 um full rescue target
        res_etd = mi.simulate_modifier_regimen("monotherapy_etd001", genotype_slc26a9="CT")
        res_etx = mi.simulate_modifier_regimen("monotherapy_etx001", genotype_slc26a9="CT")
        res_base = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")

        self.assertGreater(res_etd.hydration.steady_state_asl_height_um, res_base.hydration.steady_state_asl_height_um)
        self.assertGreater(res_etx.hydration.steady_state_asl_height_um, res_base.hydration.steady_state_asl_height_um)
        self.assertLess(res_etd.hydration.steady_state_asl_height_um, mi.PHYSIOLOGICAL_ASL_TARGET_UM)
        self.assertLess(res_etx.hydration.steady_state_asl_height_um, mi.PHYSIOLOGICAL_ASL_TARGET_UM)


class TestEpistaticSynergyAndGenotypeModifiers(unittest.TestCase):
    """Verifies SLC26A9 polymorphism epistasis and non-linear pharmacological synergy."""

    def test_slc26a9_genotype_epistasis(self):
        res_cc = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CC")
        res_ct = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        res_tt = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="TT")

        # Protective CC allele has higher baseline SLC26A9 conductance than TT
        self.assertGreater(res_cc.conductance.g_slc26a9_ms_cm2, res_ct.conductance.g_slc26a9_ms_cm2)
        self.assertGreater(res_ct.conductance.g_slc26a9_ms_cm2, res_tt.conductance.g_slc26a9_ms_cm2)

        # Protective CC allele confers higher baseline ASL height than risk TT allele
        self.assertGreater(res_cc.hydration.steady_state_asl_height_um, res_tt.hydration.steady_state_asl_height_um)

    def test_combination_epistatic_synergy(self):
        # Dual combination synergy ratio > 1.0
        res_dual = mi.simulate_modifier_regimen("dual_etd001_etx001", genotype_slc26a9="CT")
        self.assertGreater(res_dual.epistatic_synergy_ratio, 1.0)


class TestLuminalPHRecovery(unittest.TestCase):
    """Verifies bicarbonate-dependent airway surface liquid pH recovery."""

    def test_untreated_cf_acidic_asl(self):
        res_base = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        self.assertLess(res_base.hydration.luminal_ph, 6.70)

    def test_dual_combination_ph_alkalinization(self):
        res_base = mi.simulate_modifier_regimen("baseline_cf", genotype_slc26a9="CT")
        res_dual = mi.simulate_modifier_regimen("dual_etd001_etx001", genotype_slc26a9="CT")

        self.assertGreater(res_dual.hydration.luminal_ph, res_base.hydration.luminal_ph)
        self.assertGreaterEqual(res_dual.hydration.luminal_ph, mi.PHYSIOLOGICAL_PH_TARGET)


class TestAntiHallucinationAndProvenance(unittest.TestCase):
    """Verifies strict anti-hallucination provenance and prohibited clinical overclaim invariants."""

    def test_grounded_dois_and_catalog(self):
        src_text = (HERE / "modifier_interactome.py").read_text(encoding="utf-8")
        self.assertIn("10.1093/hmg/ddw290", src_text)
        self.assertIn("10.1016/j.jcf.2024.06.002", src_text)
        self.assertIn("landscape/alternative-conductance.md", src_text)

    def test_prohibited_phrases_absent(self):
        src_lower = (HERE / "modifier_interactome.py").read_text(encoding="utf-8").lower()
        prohibited = ["proven to cure", "guarantees clinical benefit", "recommended treatment regimen", "safe and effective for patients"]
        for p in prohibited:
            self.assertNotIn(p, src_lower)


class TestCLIExecution(unittest.TestCase):
    """Verifies command-line execution, self-test, and JSON export."""

    def test_cli_self_test(self):
        cmd = [sys.executable, str(HERE / "modifier_interactome.py"), "--self-test"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("PASS: modifier_interactome self-test passed.", proc.stdout)

    def test_cli_json_export(self):
        export_path = HERE / "scratch_test_modifier_export.json"
        try:
            cmd = [
                sys.executable,
                str(HERE / "modifier_interactome.py"),
                "--genotype",
                "CC",
                "--export",
                str(export_path),
                "--json",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.assertEqual(proc.returncode, 0)
            self.assertTrue(export_path.exists())

            data = json.loads(export_path.read_text(encoding="utf-8"))
            self.assertEqual(data["slc26a9_genotype"], "CC")
            self.assertIn("dual_etd001_etx001", data["regimens"])
            self.assertEqual(data["regimens"]["dual_etd001_etx001"]["verdict"], "CURATIVE_REHYDRATION_SUCCESS")
        finally:
            if export_path.exists():
                export_path.unlink()


if __name__ == "__main__":
    unittest.main()
