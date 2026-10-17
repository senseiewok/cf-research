#!/usr/bin/env python3
"""Unit tests for transcytosis_model.py (T-0177)."""
import unittest
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import transcytosis_model as tm


class TestTranscytosisModel(unittest.TestCase):
    def test_receptor_profiles_exist(self):
        self.assertIn("fcrn", tm.RECEPTOR_PROFILES)
        self.assertIn("tfr", tm.RECEPTOR_PROFILES)
        self.assertIn("pigr", tm.RECEPTOR_PROFILES)
        for rec, data in tm.RECEPTOR_PROFILES.items():
            self.assertTrue(data["pdb_id"])
            self.assertGreater(data["kd_apical_nm"], 0)
            self.assertIsInstance(data["optimal_kd_window_nm"], tuple)

    def test_fractional_occupancy(self):
        self.assertAlmostEqual(tm.calculate_fractional_occupancy(100.0, 100.0), 0.5)
        self.assertAlmostEqual(tm.calculate_fractional_occupancy(0.0, 100.0), 0.0)
        self.assertGreater(tm.calculate_fractional_occupancy(900.0, 100.0), 0.8)
        with self.assertRaises(ValueError):
            tm.calculate_fractional_occupancy(-10.0, 100.0)
        with self.assertRaises(ValueError):
            tm.calculate_fractional_occupancy(10.0, 0.0)

    def test_simulate_fcrn_optimal(self):
        res = tm.simulate_transcytosis_kinetics("fcrn")
        self.assertEqual(res["receptor_id"], "fcrn")
        self.assertEqual(res["pdb_id"], "4N0F")
        self.assertTrue(res["in_optimal_affinity_window"])
        self.assertGreaterEqual(res["transcytosis_efficiency_pct"], 40.0)
        self.assertEqual(res["verdict"], "OPTIMAL_TRANSCYTOSIS")

    def test_simulate_tfr_high_affinity_trapping(self):
        # Very high affinity (Kd = 2.0 nM) should cause lysosomal trapping
        res = tm.simulate_transcytosis_kinetics("tfr", custom_kd_nm=2.0)
        self.assertFalse(res["in_optimal_affinity_window"])
        self.assertEqual(res["verdict"], "HIGH_AFFINITY_LYSOSOMAL_TRAPPING")
        self.assertGreater(res["lysosomal_degradation_pct"], 0.0)

    def test_invalid_receptor(self):
        with self.assertRaises(ValueError):
            tm.simulate_transcytosis_kinetics("nonexistent_receptor")


if __name__ == "__main__":
    unittest.main()
