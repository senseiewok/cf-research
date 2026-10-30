#!/usr/bin/env python3
"""Comprehensive unit test suite for macrophage_evasion.py (T-0544).

Acceptance verification command:
  python -m unittest tools.discovery.test_macrophage_evasion -v && python tools/claims/check_anti_hallucination.py
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from macrophage_evasion import (
    NLM_REGISTRY_DISCLAIMER,
    PARAMETER_PROVENANCE,
    RESEARCH_BOUNDARY,
    AirwayClearanceProfile,
    AirwayVectorSimulationResult,
    CD47BindingProfile,
    calculate_cd47_binding_affinity,
    calculate_evasion_half_life,
    calculate_phagocytic_clearance_rate,
    evaluate_discovery_panel,
    get_airway_clearance_profile,
    get_cd47_binding_profile,
    self_test,
    simulate_airway_vector_kinetics,
)


class TestCD47BindingKinetics(unittest.TestCase):
    """Test Langmuir binding equilibrium for SIRP-alpha / CD47 receptor-ligand interaction."""

    def test_zero_concentration_occupancy(self):
        theta = calculate_cd47_binding_affinity(0.0, kd_uM=1.0)
        self.assertEqual(theta, 0.0)

    def test_kd_concentration_occupancy_is_fifty_percent(self):
        theta = calculate_cd47_binding_affinity(1.0, kd_uM=1.0)
        self.assertAlmostEqual(theta, 0.50, places=4)

        # Also at Kd = 0.8 uM
        theta_08 = calculate_cd47_binding_affinity(0.8, kd_uM=0.8)
        self.assertAlmostEqual(theta_08, 0.50, places=4)

    def test_high_concentration_saturation(self):
        theta = calculate_cd47_binding_affinity(100.0, kd_uM=1.0)
        self.assertGreater(theta, 0.99)
        self.assertLess(theta, 1.0)

    def test_kd_range_bounds(self):
        # K_D between 0.5 and 2.0 uM
        for kd in [0.5, 1.0, 1.5, 2.0]:
            theta = calculate_cd47_binding_affinity(4.0, kd_uM=kd)
            self.assertGreater(theta, 0.65)
            self.assertLess(theta, 0.90)

    def test_invalid_parameters_raise_value_error(self):
        with self.assertRaises(ValueError):
            calculate_cd47_binding_affinity(-1.0, kd_uM=1.0)
        with self.assertRaises(ValueError):
            calculate_cd47_binding_affinity(1.0, kd_uM=0.0)
        with self.assertRaises(ValueError):
            calculate_cd47_binding_affinity(1.0, kd_uM=-0.5)

    def test_cd47_binding_profile_dataclass(self):
        prof = get_cd47_binding_profile(4.0, kd_uM=1.0)
        self.assertIsInstance(prof, CD47BindingProfile)
        self.assertEqual(prof.cd47_conc_uM, 4.0)
        self.assertEqual(prof.kd_uM, 1.0)
        self.assertAlmostEqual(prof.receptor_occupancy_theta, 0.80, places=3)
        self.assertFalse(prof.is_saturated)

        prof_sat = get_cd47_binding_profile(20.0, kd_uM=1.0)
        self.assertTrue(prof_sat.is_saturated)


class TestPhagocyticClearanceODE(unittest.TestCase):
    """Test clearance ODE rate constants and mucosal residence half-life."""

    def test_baseline_clearance_without_cd47(self):
        # With theta = 0, clearance is k_uptake + k_ciliary
        k_clear = calculate_phagocytic_clearance_rate(k_uptake=0.35, theta_cd47=0.0, eta_evade=0.85, k_ciliary=0.05)
        self.assertAlmostEqual(k_clear, 0.40, places=4)

    def test_clearance_inhibition_with_evasion_efficiency(self):
        # With theta = 0.80 and eta_evade = 0.85
        # k_uptake_eff = 0.35 * (1 - 0.85 * 0.80) = 0.35 * (1 - 0.68) = 0.35 * 0.32 = 0.112
        # k_clear = 0.112 + 0.05 = 0.162
        k_clear = calculate_phagocytic_clearance_rate(k_uptake=0.35, theta_cd47=0.80, eta_evade=0.85, k_ciliary=0.05)
        self.assertAlmostEqual(k_clear, 0.162, places=3)

    def test_half_life_calculation(self):
        # t_1/2 = ln(2) / 0.162 approx 4.278 h
        t_half = calculate_evasion_half_life(0.162)
        self.assertAlmostEqual(t_half, math.log(2.0) / 0.162, places=3)

        # Zero clearance gives infinite half-life
        self.assertEqual(calculate_evasion_half_life(0.0), float("inf"))

    def test_half_life_threshold_benchmark(self):
        # With [CD47] = 5.0 uM, Kd = 0.8 uM, eta_evade = 0.88, k_ciliary = 0.04
        # theta = 5.0 / 5.8 = 0.8621
        # eta * theta = 0.88 * 0.8621 = 0.7586
        # k_uptake_eff = 0.35 * (1 - 0.7586) = 0.08447
        # k_clear = 0.08447 + 0.04 = 0.12447
        # t_1/2 = ln(2) / 0.12447 = 5.568 h >= 4.5 h
        prof = get_airway_clearance_profile(
            k_uptake_baseline_hr=0.35,
            cd47_conc_uM=5.0,
            kd_uM=0.8,
            eta_evade=0.88,
            k_ciliary_hr=0.04,
            target_half_life_threshold_hours=4.5,
        )
        self.assertTrue(prof.target_half_life_met)
        self.assertGreaterEqual(prof.half_life_hours, 4.5)

    def test_invalid_clearance_parameters_raise(self):
        with self.assertRaises(ValueError):
            calculate_phagocytic_clearance_rate(-0.1, theta_cd47=0.5)
        with self.assertRaises(ValueError):
            calculate_phagocytic_clearance_rate(0.35, theta_cd47=1.5)
        with self.assertRaises(ValueError):
            calculate_phagocytic_clearance_rate(0.35, theta_cd47=0.5, eta_evade=-0.1)
        with self.assertRaises(ValueError):
            calculate_phagocytic_clearance_rate(0.35, theta_cd47=0.5, k_ciliary=-0.01)


class TestAirwayVectorSimulation(unittest.TestCase):
    """Test full multi-compartment kinetic ODE simulation and mass conservation."""

    def test_mass_conservation_at_all_timepoints(self):
        res = simulate_airway_vector_kinetics(
            v_initial=100.0,
            k_uptake_baseline_hr=0.35,
            cd47_conc_uM=4.0,
            kd_uM=1.0,
            eta_evade=0.85,
            k_ciliary_hr=0.05,
            k_epithelial_hr=0.15,
            total_time_hours=24.0,
            dt_hours=0.5,
        )
        self.assertIsInstance(res, AirwayVectorSimulationResult)
        for i in range(len(res.time_points_hours)):
            total = (
                res.v_airway_fraction[i]
                + res.v_epithelial_fraction[i]
                + res.v_macrophage_fraction[i]
                + res.v_ciliary_fraction[i]
            )
            self.assertAlmostEqual(total, 100.0, delta=0.05)

    def test_monotonic_compartment_behavior(self):
        res = simulate_airway_vector_kinetics(v_initial=100.0)
        # Airway vector strictly decreases
        for i in range(1, len(res.v_airway_fraction)):
            self.assertLessEqual(res.v_airway_fraction[i], res.v_airway_fraction[i - 1])
        # Sinks strictly increase
        for i in range(1, len(res.v_epithelial_fraction)):
            self.assertGreaterEqual(res.v_epithelial_fraction[i], res.v_epithelial_fraction[i - 1])
            self.assertGreaterEqual(res.v_macrophage_fraction[i], res.v_macrophage_fraction[i - 1])
            self.assertGreaterEqual(res.v_ciliary_fraction[i], res.v_ciliary_fraction[i - 1])

    def test_fold_enhancement_over_bare_vector(self):
        res_stealth = simulate_airway_vector_kinetics(cd47_conc_uM=5.0, kd_uM=0.8, eta_evade=0.88)
        self.assertGreaterEqual(res_stealth.fold_enhancement_vs_bare, 1.50)
        self.assertTrue(res_stealth.curative_delivery_benchmark_met)


class TestDiscoveryPanelAndProvenance(unittest.TestCase):
    """Test panel evaluation, AH-9 parameter provenance, and regulatory compliance."""

    def test_panel_evaluation_structure(self):
        panel = evaluate_discovery_panel()
        self.assertEqual(panel["status"], "PASS")
        self.assertEqual(panel["trial_surveillance"], "NCT05668741")
        self.assertEqual(panel["nlm_attribution"], NLM_REGISTRY_DISCLAIMER)
        self.assertEqual(panel["research_boundary"], RESEARCH_BOUNDARY)
        self.assertEqual(panel["configurations_evaluated"], 4)

        # Increasing CD47 density monotonically increases half-life
        configs = panel["panel_results"]
        for i in range(1, len(configs)):
            self.assertGreater(configs[i]["half_life_hours"], configs[i - 1]["half_life_hours"])

    def test_self_test_function(self):
        self.assertTrue(self_test())

    def test_parameter_provenance_bounds(self):
        for name, meta in PARAMETER_PROVENANCE.items():
            val = meta["value"]
            low, high = meta["bounds"]
            self.assertGreaterEqual(val, low)
            self.assertLessEqual(val, high)
            self.assertIn(meta["source_type"], ["empirical", "derived", "phenomenological_heuristic"])

    def test_cli_execution_with_json(self):
        script = HERE / "macrophage_evasion.py"
        cmd = [sys.executable, str(script), "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        data = json.loads(proc.stdout)
        self.assertEqual(data["status"], "PASS")
        self.assertEqual(data["trial_surveillance"], "NCT05668741")


if __name__ == "__main__":
    unittest.main()
