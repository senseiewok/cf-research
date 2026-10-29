#!/usr/bin/env python3
"""Unit tests for Inhaled mRNA-LNP Endosomal Escape Kinetics & Buffer Capacity Model (T-0540)."""
from __future__ import annotations

import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from tools.discovery.lnp_endosomal_escape import (
    FORMULATION_ARCHETYPES,
    LNP_ESCAPE_PARAMETERS,
    calculate_ionizable_lipid_protonation,
    calculate_volumetric_buffer_capacity,
    evaluate_falsification_negative_control,
    self_test,
    simulate_endosomal_acidification_and_escape,
)


class TestIonizableLipidProtonation(unittest.TestCase):
    """Test Henderson-Hasselbalch protonation calculations."""

    def test_ph_equals_pka(self) -> None:
        alpha = calculate_ionizable_lipid_protonation(ph=6.45, pka=6.45)
        self.assertAlmostEqual(alpha, 0.50, places=3)

    def test_acidic_ph_full_protonation(self) -> None:
        alpha = calculate_ionizable_lipid_protonation(ph=4.5, pka=6.45)
        self.assertGreater(alpha, 0.95)

    def test_basic_ph_deprotonated(self) -> None:
        alpha = calculate_ionizable_lipid_protonation(ph=8.5, pka=6.45)
        self.assertLess(alpha, 0.05)

    def test_invalid_ph_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_ionizable_lipid_protonation(ph=-1.0, pka=6.45)
        with self.assertRaises(ValueError):
            calculate_ionizable_lipid_protonation(ph=15.0, pka=6.45)

    def test_invalid_pka_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_ionizable_lipid_protonation(ph=7.0, pka=-0.5)


class TestBufferCapacity(unittest.TestCase):
    """Test volumetric buffer capacity (proton sponge effect)."""

    def test_peak_at_pka(self) -> None:
        pka = 6.45
        conc = 25.0
        beta_peak = calculate_volumetric_buffer_capacity(ph=pka, pka=pka, total_lipid_conc_mm=conc)
        beta_low = calculate_volumetric_buffer_capacity(ph=5.0, pka=pka, total_lipid_conc_mm=conc)
        beta_high = calculate_volumetric_buffer_capacity(ph=8.0, pka=pka, total_lipid_conc_mm=conc)
        self.assertGreater(beta_peak, beta_low)
        self.assertGreater(beta_peak, beta_high)
        # Theoretical max = 2.302585 * 25.0 / 4 = 14.391 mM/pH
        self.assertAlmostEqual(beta_peak, 2.302585 * 25.0 / 4.0, places=2)

    def test_invalid_lipid_conc_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            calculate_volumetric_buffer_capacity(ph=6.45, pka=6.45, total_lipid_conc_mm=0.0)


class TestEndosomalEscapeSimulation(unittest.TestCase):
    """Test ODE simulation across formulations."""

    def test_vx522_archetype_meets_therapeutic_threshold(self) -> None:
        res = simulate_endosomal_acidification_and_escape("vx522_archetype")
        self.assertEqual(res.formulation_id, "vx522_archetype")
        self.assertTrue(res.therapeutic_threshold_met)
        self.assertGreaterEqual(res.cytosolic_escape_efficiency_pct, 1.5)
        self.assertLessEqual(res.cytosolic_escape_efficiency_pct, 3.5)
        self.assertTrue(res.hexagonal_phase_formed)
        self.assertGreater(res.cytosolic_mrna_released_pct, 50.0)

    def test_rct2100_archetype_meets_therapeutic_threshold(self) -> None:
        res = simulate_endosomal_acidification_and_escape("rct2100_archetype")
        self.assertEqual(res.formulation_id, "rct2100_archetype")
        self.assertTrue(res.therapeutic_threshold_met)
        self.assertGreaterEqual(res.cytosolic_escape_efficiency_pct, 1.5)
        self.assertLessEqual(res.cytosolic_escape_efficiency_pct, 3.5)
        self.assertTrue(res.hexagonal_phase_formed)

    def test_neutral_control_fails_therapeutic_threshold(self) -> None:
        res = simulate_endosomal_acidification_and_escape("neutral_non_ionizable_control")
        self.assertFalse(res.therapeutic_threshold_met)
        self.assertLess(res.cytosolic_escape_efficiency_pct, 0.50)
        self.assertFalse(res.hexagonal_phase_formed)

    def test_unknown_formulation_raises_key_error(self) -> None:
        with self.assertRaises(KeyError):
            simulate_endosomal_acidification_and_escape("unknown_lnp")

    def test_invalid_transit_time_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            simulate_endosomal_acidification_and_escape("vx522_archetype", transit_duration_min=0.0)


class TestProvenanceAndFalsification(unittest.TestCase):
    """Test AH-9 parameter provenance and AH-10 negative controls."""

    def test_ah9_parameter_provenance(self) -> None:
        self.assertGreaterEqual(len(LNP_ESCAPE_PARAMETERS), 3)
        for p_name, param in LNP_ESCAPE_PARAMETERS.items():
            self.assertTrue(len(param.doi_or_pdb) > 0, f"Missing DOI for {p_name}")
            self.assertIn(param.category, ("empirical", "derived", "phenomenological_heuristic"))

    def test_ah10_falsification_negative_control(self) -> None:
        neg = evaluate_falsification_negative_control()
        self.assertEqual(neg["rule"], "AH-10")
        self.assertTrue(neg["passed"])

    def test_self_test_function(self) -> None:
        ret = self_test()
        self.assertEqual(ret, 0)

    def test_cli_execution_and_json(self) -> None:
        cmd = [sys.executable, str(HERE / "lnp_endosomal_escape.py"), "--formulation", "vx522_archetype", "--json"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        data = json.loads(proc.stdout)
        self.assertEqual(data["formulation_id"], "vx522_archetype")
        self.assertTrue(data["therapeutic_threshold_met"])


if __name__ == "__main__":
    unittest.main()
