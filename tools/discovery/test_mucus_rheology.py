#!/usr/bin/env python3
"""Deterministic verification test suite for mucus_rheology.py (T-0185).

Dual compatibility:
- pytest tools/discovery/test_mucus_rheology.py
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

import mucus_rheology as mr


class TestViscoelasticModulusRheology(unittest.TestCase):
    """Verifies complex shear modulus G*(omega), G', G'', and loss tangent scaling."""

    def test_g_prime_monotonic_increase(self):
        # Storage modulus should strictly increase with solids content
        solids = [1.5, 2.0, 3.5, 5.0, 7.5, 10.0]
        g_primes = [mr.calculate_viscoelastic_modulus(s).storage_modulus_g_prime_pa for s in solids]
        for i in range(len(g_primes) - 1):
            self.assertLess(g_primes[i], g_primes[i + 1])

    def test_loss_tangent_transition(self):
        # Healthy 2% mucus is dissipative sol/weak gel (tan delta ~ 0.8-1.2)
        mod_healthy = mr.calculate_viscoelastic_modulus(2.0)
        self.assertGreater(mod_healthy.loss_tangent_tan_delta, 0.70)
        self.assertIn(mod_healthy.gel_state_regime, ("DISSIPATIVE_SOL_FLUID", "WEAK_VISCOELASTIC_GEL"))

        # Severe CF 8% mucus is a stiff elastic crosslinked gel (tan delta < 0.35)
        mod_cf = mr.calculate_viscoelastic_modulus(8.0)
        self.assertLess(mod_cf.loss_tangent_tan_delta, 0.35)
        self.assertEqual(mod_cf.gel_state_regime, "ELASTIC_CROSS_LINKED_GEL")

    def test_zero_shear_viscosity_scaling(self):
        mod_healthy = mr.calculate_viscoelastic_modulus(2.0)
        mod_cf = mr.calculate_viscoelastic_modulus(8.0)
        # CF zero-shear viscosity should exceed healthy by at least 20-fold
        self.assertGreater(mod_cf.zero_shear_viscosity_pa_s, mod_healthy.zero_shear_viscosity_pa_s * 20.0)

    def test_frequency_dependence(self):
        # Higher oscillatory shear frequency increases complex modulus magnitude
        mod_low_freq = mr.calculate_viscoelastic_modulus(4.0, frequency_rad_s=0.5)
        mod_high_freq = mr.calculate_viscoelastic_modulus(4.0, frequency_rad_s=5.0)
        self.assertGreater(mod_high_freq.complex_modulus_magnitude_pa, mod_low_freq.complex_modulus_magnitude_pa)

    def test_invalid_solids_or_frequency(self):
        with self.assertRaises(ValueError):
            mr.calculate_viscoelastic_modulus(0.2)  # Below physiological range
        with self.assertRaises(ValueError):
            mr.calculate_viscoelastic_modulus(25.0)  # Above physiological range
        with self.assertRaises(ValueError):
            mr.calculate_viscoelastic_modulus(2.0, frequency_rad_s=-1.0)


class TestMucinPoreScalingAndObstruction(unittest.TestCase):
    """Verifies MUC5B/MUC5AC pore dimension scaling xi(w) and steric boundaries."""

    def test_pore_size_scaling_bounds(self):
        pore_healthy = mr.calculate_mucin_pore_size(2.0)
        pore_cf = mr.calculate_mucin_pore_size(8.0)

        # Healthy pore size >= 160 nm
        self.assertGreaterEqual(pore_healthy.pore_diameter_nm, 160.0)
        self.assertEqual(pore_healthy.regime_description, "EXPANDED_HEALTHY_MESH_PERMISSIVE")

        # CF pore size collapses below 75 nm
        self.assertLessEqual(pore_cf.pore_diameter_nm, 75.0)
        self.assertEqual(pore_cf.regime_description, "COLLAPSED_CF_MESH_STERICALLY_OCCLUDED")

    def test_mucin_volume_fraction_scaling(self):
        pore_healthy = mr.calculate_mucin_pore_size(2.0)
        pore_cf = mr.calculate_mucin_pore_size(8.0)
        self.assertLess(pore_healthy.mucin_volume_fraction, pore_cf.mucin_volume_fraction)


class TestHinderedDiffusivity(unittest.TestCase):
    """Verifies Ogston-Mackie steric obstruction and surface chemistry adhesion penalties."""

    def test_peg_vs_cationic_diffusivity(self):
        # At 40 nm in 6% mucus, PEG-brush diffuses significantly faster than unshielded cationic
        res_peg = mr.calculate_hindered_particle_diffusivity(40.0, "peg_brush", solids_pct=6.0)
        res_cat = mr.calculate_hindered_particle_diffusivity(40.0, "unshielded_cationic", solids_pct=6.0)

        self.assertGreater(res_peg.relative_diffusivity_ratio, res_cat.relative_diffusivity_ratio * 30.0)
        self.assertEqual(res_peg.surface_adhesion_penalty, 1.0)
        self.assertLess(res_cat.surface_adhesion_penalty, 0.05)

    def test_steric_mesh_cutoff_entrapment(self):
        # 100 nm particle in 8% mucus (pore ~55 nm) undergoes severe steric entrapment
        res_large = mr.calculate_hindered_particle_diffusivity(100.0, "peg_brush", solids_pct=8.0)
        self.assertLessEqual(res_large.relative_diffusivity_ratio, 0.001)
        self.assertEqual(res_large.transmucosal_penetration_status, "SEVERE_MUCIN_ENTRAPMENT")

    def test_small_aso_rapid_percolation(self):
        # Naked ASO (~2.5 nm) rapidly permeates even 8% CF mucus
        res_aso = mr.calculate_hindered_particle_diffusivity(2.5, "aso_oligo", solids_pct=8.0)
        self.assertGreater(res_aso.relative_diffusivity_ratio, 0.45)
        self.assertEqual(res_aso.transmucosal_penetration_status, "RAPID_TRANSMUCOSAL_PERMEATION")


class TestMucociliaryTransportAndRescue(unittest.TestCase):
    """Verifies mucociliary transport velocity v_MCT and mucolytic rescue interventions."""

    def test_transport_velocity_collapse_in_cf(self):
        mct_healthy = mr.calculate_mucociliary_transport_velocity(2.0)
        mct_cf = mr.calculate_mucociliary_transport_velocity(8.0)

        # Normal clearance velocity > 50 um/s
        self.assertGreater(mct_healthy.transport_velocity_um_s, 50.0)
        self.assertEqual(mct_healthy.therapeutic_clearance_competence, "PHYSIOLOGICAL_HIGH_CLEARANCE")

        # CF clearance collapses below 10 um/s (mucostasis)
        self.assertLess(mct_cf.transport_velocity_um_s, 10.0)
        self.assertEqual(mct_cf.therapeutic_clearance_competence, "MUCOSTASIS_CLEARANCE_FAILURE")

    def test_mucolytic_rescue_delta(self):
        # Hypertonic saline 7% rehydrates CF mucus, boosting v_MCT significantly
        rescue_hs = mr.simulate_mucolytic_rescue(8.0, "hypertonic_saline_7pct")
        self.assertGreater(rescue_hs.transport_fold_improvement, 3.0)
        self.assertLess(rescue_hs.post_treatment_solids_pct, rescue_hs.initial_solids_pct)
        self.assertIn(rescue_hs.verdict, ("PARTIAL_FUNCTIONAL_CLEARANCE_RESTORED", "COMPLETE_PHYSIOLOGICAL_CLEARANCE_RESTORED"))

    def test_triple_combination_synergy(self):
        # Triple combination achieves greater improvement than monotherapy
        rescue_hs = mr.simulate_mucolytic_rescue(8.0, "hypertonic_saline_7pct")
        rescue_triple = mr.simulate_mucolytic_rescue(8.0, "triple_combination")
        self.assertGreater(rescue_triple.transport_fold_improvement, rescue_hs.transport_fold_improvement)
        self.assertEqual(rescue_triple.verdict, "COMPLETE_PHYSIOLOGICAL_CLEARANCE_RESTORED")


class TestAntiHallucinationAndProvenance(unittest.TestCase):
    """Verifies strict anti-hallucination provenance and prohibited clinical overclaim invariants."""

    def test_nct_provenance_in_snapshot(self):
        snapshot_file = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
        self.assertTrue(snapshot_file.exists())
        with open(snapshot_file, encoding="utf-8") as f:
            snap_data = json.load(f)
        valid_ncts = {s["nct_id"] for s in snap_data.get("studies", []) if "nct_id" in s}

        # Assert NCT05668741 (VX-522) is present in snapshot
        self.assertIn("NCT05668741", valid_ncts)

    def test_prohibited_phrases_absent(self):
        src_text = (HERE / "mucus_rheology.py").read_text(encoding="utf-8").lower()
        prohibited = ["proven to cure", "guarantees clinical benefit", "recommended treatment regimen", "safe and effective for patients"]
        for p in prohibited:
            self.assertNotIn(p, src_text)


class TestCLIExecution(unittest.TestCase):
    """Verifies command-line execution, self-test, and JSON export."""

    def test_cli_self_test(self):
        cmd = [sys.executable, str(HERE / "mucus_rheology.py"), "--self-test"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("PASS: mucus_rheology self-test passed.", proc.stdout)

    def test_cli_json_export(self):
        export_path = HERE / "scratch_test_rheology_export.json"
        try:
            cmd = [
                sys.executable,
                str(HERE / "mucus_rheology.py"),
                "--solids-pct",
                "7.0",
                "--export",
                str(export_path),
                "--json",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.assertEqual(proc.returncode, 0)
            self.assertTrue(export_path.exists())

            data = json.loads(export_path.read_text(encoding="utf-8"))
            self.assertEqual(data["provenance"]["nct_id"], "NCT05668741")
            self.assertIn("viscoelastic_rheology", data)
            self.assertIn("mucin_pore_scaling", data)
            self.assertIn("particle_diffusivity", data)
            self.assertIn("mucociliary_transport", data)
        finally:
            if export_path.exists():
                export_path.unlink()


if __name__ == "__main__":
    unittest.main()
