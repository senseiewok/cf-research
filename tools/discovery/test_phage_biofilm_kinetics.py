#!/usr/bin/env python3
"""Deterministic verification test suite for phage_biofilm_kinetics.py (T-0197).

Dual test runner compatibility:
- pytest tools/discovery/test_phage_biofilm_kinetics.py -v
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

import phage_biofilm_kinetics as pbk


class TestLyticReplicationKinetics(unittest.TestCase):
    """Verifies differential delay replication dynamics and lytic decimation."""

    def setUp(self):
        self.cocktail = pbk.get_standard_phage_cocktail()
        self.omko1 = self.cocktail["Phage_OMKO1_like"]

    def test_lytic_killing_monotonic_reduction(self):
        snapshots = pbk.calculate_lytic_kinetics(
            initial_bacteria_cfu_ml=1.0e8,
            initial_phage_pfu_ml=1.0e9,
            phage=self.omko1,
            duration_hours=6.0,
            time_step_hours=0.1,
        )
        self.assertGreater(len(snapshots), 50)
        # End bacterial load must be substantially lower than initial load
        final = snapshots[-1]
        self.assertGreater(final.log10_bacterial_reduction, 2.0)
        self.assertGreater(final.lysed_cells_cumulative_ml, 1.0e7)

    def test_pas_sub_mic_accelerates_clearance(self):
        # Simulation without antibiotic
        sim_no_abx = pbk.calculate_lytic_kinetics(
            initial_bacteria_cfu_ml=1.0e8,
            initial_phage_pfu_ml=5.0e8,
            phage=self.omko1,
            duration_hours=4.0,
            sub_mic_abx_ratio=0.0,
        )
        # Simulation with sub-MIC antibiotic (PAS benefit)
        sim_pas = pbk.calculate_lytic_kinetics(
            initial_bacteria_cfu_ml=1.0e8,
            initial_phage_pfu_ml=5.0e8,
            phage=self.omko1,
            duration_hours=4.0,
            sub_mic_abx_ratio=0.35,
        )
        # PAS-enhanced replication achieves higher log10 bacterial reduction
        self.assertGreaterEqual(
            sim_pas[-1].log10_bacterial_reduction,
            sim_no_abx[-1].log10_bacterial_reduction,
        )

    def test_invalid_input_parameters(self):
        with self.assertRaises(ValueError):
            pbk.calculate_lytic_kinetics(
                initial_bacteria_cfu_ml=-100,
                initial_phage_pfu_ml=1.0e8,
                phage=self.omko1,
            )
        with self.assertRaises(ValueError):
            pbk.calculate_lytic_kinetics(
                initial_bacteria_cfu_ml=1.0e8,
                initial_phage_pfu_ml=0,
                phage=self.omko1,
            )
        with self.assertRaises(ValueError):
            pbk.calculate_lytic_kinetics(
                initial_bacteria_cfu_ml=1.0e8,
                initial_phage_pfu_ml=1.0e8,
                phage=self.omko1,
                duration_hours=-1.0,
            )


class TestBiofilmHindranceAndLyase(unittest.TestCase):
    """Verifies alginate degradation and restored hindered diffusivity in biofilm."""

    def setUp(self):
        self.cocktail = pbk.get_standard_phage_cocktail()
        self.lyase_phage = self.cocktail["Phage_PaWRA_lyase"]
        self.non_lyase_phage = self.cocktail["Phage_OMKO1_like"]

    def test_alginate_steric_hindrance_without_lyase(self):
        rem_pct, d_eff, hindrance = pbk.calculate_biofilm_hindrance(
            alginate_density_mg_ml=7.5,
            has_lyase=False,
            lyase_activity_u_ml=0.0,
            elapsed_minutes=45.0,
        )
        self.assertEqual(rem_pct, 100.0)
        self.assertEqual(hindrance, pbk.UNHYDROLYZED_BIOFILM_HINDRANCE_FACTOR)
        self.assertLess(d_eff, pbk.DIFFUSIVITY_WATER_PHAGE_CM2_S * 0.05)

    def test_alginate_depolymerization_with_lyase(self):
        rem_pct, d_eff, hindrance = pbk.calculate_biofilm_hindrance(
            alginate_density_mg_ml=7.5,
            has_lyase=True,
            lyase_activity_u_ml=25.0,
            elapsed_minutes=60.0,
        )
        # Alginate should be significantly hydrolyzed
        self.assertLess(rem_pct, 50.0)
        # Diffusivity should be restored by several-fold
        self.assertGreater(hindrance, pbk.UNHYDROLYZED_BIOFILM_HINDRANCE_FACTOR * 5.0)
        self.assertGreater(d_eff, pbk.DIFFUSIVITY_WATER_PHAGE_CM2_S * 0.20)

    def test_depth_penetration_gradient(self):
        profiles = pbk.simulate_biofilm_depth_penetration(
            biofilm_thickness_um=50.0,
            initial_phage_titer_pfu_ml=1.0e9,
            phage=self.lyase_phage,
            time_minutes=45.0,
            depth_slices=5,
        )
        self.assertEqual(len(profiles), 6)
        # Surface concentration should be highest
        self.assertEqual(profiles[0].depth_um, 0.0)
        self.assertGreater(profiles[0].local_phage_titer_pfu_ml, profiles[-1].local_phage_titer_pfu_ml)
        # Bacterial kill efficiency should decline with depth
        self.assertGreaterEqual(
            profiles[0].bacterial_kill_efficiency_pct,
            profiles[-1].bacterial_kill_efficiency_pct,
        )


class TestEvolutionarySteeringAndTradeoffs(unittest.TestCase):
    """Verifies MexAB-OprM efflux collapse, fitness costs, and resensitization."""

    def test_oprm_resensitization_and_fitness_cost(self):
        oprm = pbk.evaluate_evolutionary_steering("OprM")
        self.assertIn("OprM", oprm.target_receptor)
        # 8-fold ciprofloxacin resensitization (4.0 -> 0.5 ug/mL)
        self.assertEqual(oprm.resensitization_fold_cipro, 8.0)
        self.assertEqual(oprm.ciprofloxacin_mic_ug_ml, 0.5)
        # 8-fold ceftazidime resensitization (32.0 -> 4.0 ug/mL)
        self.assertEqual(oprm.resensitization_fold_ceftaz, 8.0)
        self.assertEqual(oprm.ceftazidime_mic_ug_ml, 4.0)
        # Fitness cost ~ 18% growth reduction
        self.assertEqual(oprm.fitness_cost_pct, 18.0)
        self.assertLess(oprm.bacterial_growth_rate_per_h, pbk.DEFAULT_BACTERIAL_GROWTH_RATE_MAX_PER_H)

    def test_pila_biofilm_cripple(self):
        pila = pbk.evaluate_evolutionary_steering("PilA")
        self.assertIn("PilA", pila.target_receptor)
        self.assertEqual(pila.fitness_cost_pct, 8.0)
        self.assertIn("twitching motility", pila.phenotypic_tradeoff)

    def test_lps_rough_transition(self):
        lps = pbk.evaluate_evolutionary_steering("LPS")
        self.assertIn("LPS", lps.target_receptor)
        self.assertEqual(lps.fitness_cost_pct, 15.0)
        self.assertIn("serum resistance", lps.phenotypic_tradeoff)

    def test_unknown_receptor_raises(self):
        with self.assertRaises(ValueError):
            pbk.evaluate_evolutionary_steering("UnknownReceptorX")


class TestPhageAntibioticSynergy(unittest.TestCase):
    """Verifies sub-MIC enhancement, burst acceleration, and FICI classification."""

    def setUp(self):
        self.cocktail = pbk.get_standard_phage_cocktail()
        self.omko1 = self.cocktail["Phage_OMKO1_like"]

    def test_sub_mic_ciprofloxacin_synergy(self):
        # 1.0 ug/mL ciprofloxacin is 0.25 * MIC (sub-MIC)
        res = pbk.calculate_pas_synergy("Ciprofloxacin", 1.0, self.omko1)
        self.assertTrue(res.is_sub_mic)
        self.assertEqual(res.fractional_mic_ratio, 0.25)
        self.assertGreater(res.burst_multiplier, 1.2)
        self.assertGreater(res.effective_burst_size_pfu, self.omko1.burst_size_pfu)
        self.assertLessEqual(res.fici_score, 0.5)
        self.assertEqual(res.synergy_classification, "SYNERGISTIC")
        self.assertGreater(res.synergy_excess_log10, 1.0)

    def test_high_dose_inhibits_phage_burst(self):
        # 8.0 ug/mL ciprofloxacin is 2.0 * MIC (supra-inhibitory)
        res = pbk.calculate_pas_synergy("Ciprofloxacin", 8.0, self.omko1)
        self.assertFalse(res.is_sub_mic)
        self.assertGreater(res.fractional_mic_ratio, 1.0)
        # Excessive antibiotic halts metabolism and suppresses burst
        self.assertLess(res.burst_multiplier, 0.5)


class TestCocktailVsMonotherapy(unittest.TestCase):
    """Verifies multi-receptor cocktail suppression of escape mutations."""

    def test_cocktail_resistance_probability_advantage(self):
        tri = pbk.simulate_cocktail_vs_monotherapy(use_tri_cocktail=True)
        mono = pbk.simulate_cocktail_vs_monotherapy(use_tri_cocktail=False)

        # Monotherapy escape probability is 10^-6
        self.assertEqual(mono.resistance_escape_probability, 1.0e-6)
        # Tri-modal cocktail escape probability is 10^-18
        self.assertAlmostEqual(tri.resistance_escape_probability, 1.0e-18)

        # Cocktail yields superior total decimation and biomass reduction
        self.assertGreater(tri.total_log10_reduction, mono.total_log10_reduction + 2.0)
        self.assertGreater(tri.biofilm_biomass_reduction_pct, 90.0)

    def test_trial_grounding_anchors(self):
        tri = pbk.simulate_cocktail_vs_monotherapy(use_tri_cocktail=True)
        self.assertIn("NCT04684641", tri.clinical_trial_grounding)
        self.assertIn("NCT05453578", tri.clinical_trial_grounding)
        self.assertIn("NCT05616221", tri.clinical_trial_grounding)


class TestManifestAndCLI(unittest.TestCase):
    """Verifies complete manifest generation, self-checks, and CLI functionality."""

    def test_manifest_completeness(self):
        manifest = pbk.export_complete_phage_biophysics_manifest()
        self.assertEqual(manifest["task_id"], "T-0197")
        self.assertEqual(manifest["answer_tier"], "Tier T2")
        self.assertEqual(len(manifest["clinical_trial_anchors"]), 3)
        self.assertEqual(len(manifest["primary_dois"]), 4)
        self.assertIn("NCT04684641", manifest["clinical_trial_anchors"])
        self.assertIn("10.1038/s41591-021-01659-5", manifest["primary_dois"])
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", manifest["registry_disclaimer"])

    def test_cli_check_flag(self):
        cmd = [sys.executable, str(HERE / "phage_biofilm_kinetics.py"), "--check"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("PASS: Self-diagnostic consistency checks verified", result.stdout)

    def test_cli_json_flag(self):
        cmd = [sys.executable, str(HERE / "phage_biofilm_kinetics.py"), "--json"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        data = json.loads(result.stdout)
        self.assertEqual(data["task_id"], "T-0197")


if __name__ == "__main__":
    unittest.main()
