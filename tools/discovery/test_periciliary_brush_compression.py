"""Unit tests for Periciliary Brush Osmotic Compression vs Airway Surface Liquid Height (T-0562).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1126/science.1223012 (A Periciliary Brush Promotes the Lung Health by Separating the Mucus Layer from Airway Epithelia)
  - DOI:10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)
  - DOI:10.1016/j.biomaterials.2008.12.076 (The penetration of fresh undiluted sputum expectorated by cystic fibrosis patients)
  - NCT05668741 (Phase 1/2 VX-522 inhaled mRNA LNP study in CF)
"""
from __future__ import annotations

import unittest
from tools.discovery.periciliary_brush_compression import (
    PericiliaryBrushParameters,
    calculate_equilibrium_pcl_height,
    compute_mucus_osmotic_pressure,
    compute_pcl_osmotic_pressure,
    evaluate_periciliary_brush_compression,
    falsify_periciliary_brush_compression,
    simulate_periciliary_brush,
)


class PericiliaryBrushCompressionTests(unittest.TestCase):
    """Test suite verifying Flory-Huggins grafted brush mechanics and compressive collapse."""

    def setUp(self) -> None:
        self.report = evaluate_periciliary_brush_compression()

    def test_report_structure_and_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0562")
        self.assertTrue(self.report["periciliary_brush_verified"])
        self.assertIn("Healthy_Hydrated", self.report["profiles"])
        self.assertIn("CF_Untreated_Hyperconcentrated", self.report["profiles"])
        self.assertIn("CF_Partially_Rescued", self.report["profiles"])

    def test_mucus_osmotic_pressure_scaling(self) -> None:
        # Pi_mucus = A * (phi / phi_0)^B
        # At baseline phi = 2.0%, A = 200, B = 2.5: Pi_mucus = 200.0 Pa
        pi_baseline = compute_mucus_osmotic_pressure(2.0, baseline_pct=2.0, coeff_a=200.0, exp_b=2.5)
        self.assertAlmostEqual(pi_baseline, 200.0, places=2)

        # At phi = 4.0% (2x baseline): Pi_mucus = 200 * (2^2.5) = 200 * 5.65685 = 1131.37 Pa
        pi_4pct = compute_mucus_osmotic_pressure(4.0, baseline_pct=2.0, coeff_a=200.0, exp_b=2.5)
        self.assertAlmostEqual(pi_4pct, 1131.37, places=1)

        # At phi = 8.0% (4x baseline): Pi_mucus = 200 * (4^2.5) = 200 * 32 = 6400.0 Pa
        pi_8pct = compute_mucus_osmotic_pressure(8.0, baseline_pct=2.0, coeff_a=200.0, exp_b=2.5)
        self.assertAlmostEqual(pi_8pct, 6400.0, places=1)

    def test_pcl_osmotic_pressure_flory_huggins(self) -> None:
        # Pi_PCL = Pi_0 * (h_0 / h)^(9/4)
        # At h = h_0 = 7.0 um: Pi_PCL = Pi_0 = 500.0 Pa
        pi_uncompressed = compute_pcl_osmotic_pressure(7.0, resting_pi_0=500.0, uncompressed_h_0=7.0)
        self.assertAlmostEqual(pi_uncompressed, 500.0, places=2)

        # At compressed h = 3.5 um (half height): Pi_PCL = 500 * (2^2.25) = 500 * 4.7568 = 2378.41 Pa
        pi_compressed = compute_pcl_osmotic_pressure(3.5, resting_pi_0=500.0, uncompressed_h_0=7.0)
        self.assertAlmostEqual(pi_compressed, 2378.41, places=1)

    def test_equilibrium_height_and_motility_stratification(self) -> None:
        healthy = self.report["profiles"]["Healthy_Hydrated"]
        cf_untreated = self.report["profiles"]["CF_Untreated_Hyperconcentrated"]
        cf_rescued = self.report["profiles"]["CF_Partially_Rescued"]
        severe_cf = self.report["profiles"]["Severe_CF_Dehydrated"]

        # Healthy: full uncompressed height (7.0 um) and normal CBF (12.0 Hz)
        self.assertEqual(healthy["equilibrium_height_um"], 7.0)
        self.assertEqual(healthy["compression_ratio"], 1.0)
        self.assertEqual(healthy["ciliary_motility_status"], "ACTIVE_UNCOMPRESSED_MOTILITY")
        self.assertEqual(healthy["effective_cbf_hz"], 12.0)

        # CF untreated: collapsed PCL height < 4.0 um (critical collapse) and ciliary arrest
        self.assertLess(cf_untreated["equilibrium_height_um"], 4.0)
        self.assertGreater(cf_untreated["equilibrium_height_um"], 1.5)
        self.assertEqual(cf_untreated["ciliary_motility_status"], "ARRESTED_CILIARY_COLLAPSE")
        self.assertLess(cf_untreated["effective_cbf_hz"], 4.0)

        # Partially rescued: intermediate PCL height (between 4.0 and 7.0 um)
        self.assertGreaterEqual(cf_rescued["equilibrium_height_um"], 4.0)
        self.assertLess(cf_rescued["equilibrium_height_um"], 7.0)
        self.assertEqual(cf_rescued["ciliary_motility_status"], "COMPROMISED_PARTIAL_COMPRESSION")
        self.assertGreater(cf_rescued["effective_cbf_hz"], cf_untreated["effective_cbf_hz"])

        # Severe CF: extreme collapse
        self.assertLess(severe_cf["equilibrium_height_um"], cf_untreated["equilibrium_height_um"])

    def test_ah10_dynamic_falsification_negative_control(self) -> None:
        fals_res = falsify_periciliary_brush_compression()
        self.assertEqual(fals_res["rule"], "AH-10")
        self.assertEqual(fals_res["verdict"], "PASS")
        self.assertTrue(fals_res["negative_control_passed"])
        self.assertLess(fals_res["collapsed_equilibrium_height_um"], 2.0)
        self.assertEqual(fals_res["effective_cbf_hz"], 0.0)


if __name__ == "__main__":
    unittest.main()
