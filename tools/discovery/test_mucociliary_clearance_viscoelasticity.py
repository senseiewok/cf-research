"""Unit tests for Mucociliary Transport Velocity vs Viscoelastic Modulus & CBF (T-0558).

Primary Biophysical Grounding & Evidence:
  - DOI:10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)
  - DOI:10.1016/j.biomaterials.2008.12.076 (The penetration of fresh undiluted sputum expectorated by cystic fibrosis patients)
  - DOI:10.1126/science.1223012 (A Periciliary Brush Promotes Lung Health by Separating the Mucus Layer)
  - NCT05668741 (Phase 1/2 VX-522 inhaled mRNA LNP study in CF)
"""
from __future__ import annotations

import unittest
from tools.discovery.mucociliary_clearance_viscoelasticity import (
    AirwayMucusRheologyParameters,
    compute_complex_shear_modulus,
    evaluate_mucociliary_clearance,
    falsify_mucociliary_clearance,
    simulate_mucociliary_transport,
)


class MucociliaryClearanceViscoelasticityTests(unittest.TestCase):
    """Test suite verifying ciliary forward stroke mechanics, mucus complex shear modulus, and clearance velocity."""

    def setUp(self) -> None:
        self.report = evaluate_mucociliary_clearance()

    def test_report_structure_and_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0558")
        self.assertTrue(self.report["mucociliary_clearance_verified"])
        self.assertIn("Healthy_Hydrated", self.report["profiles"])
        self.assertIn("CF_Untreated_Hyperconcentrated", self.report["profiles"])

    def test_complex_shear_modulus_computation(self) -> None:
        # |G*| = sqrt((G')^2 + (G'')^2)
        g_prime = 3.0
        g_double_prime = 4.0
        self.assertAlmostEqual(compute_complex_shear_modulus(g_prime, g_double_prime), 5.0, places=5)
        self.assertEqual(compute_complex_shear_modulus(0.0, 0.0), 0.0)

    def test_clearance_velocity_stratification(self) -> None:
        healthy = self.report["profiles"]["Healthy_Hydrated"]
        cf_untreated = self.report["profiles"]["CF_Untreated_Hyperconcentrated"]
        cf_rescued = self.report["profiles"]["CF_Partially_Hydrated"]
        pcd = self.report["profiles"]["Primary_Ciliary_Dyskinesia_Control"]

        # Healthy: active clearance (>= 40 um/s, equivalent to >= 2.4 mm/min)
        self.assertGreaterEqual(healthy["transport_velocity_um_per_s"], 40.0)
        self.assertEqual(healthy["clearance_status"], "ACTIVE_CLEARANCE")
        self.assertEqual(healthy["pcl_brush_compression_ratio"], 1.0)

        # CF untreated: arrested mucus stasis (< 1.0 um/s)
        self.assertLess(cf_untreated["transport_velocity_um_per_s"], 1.0)
        self.assertEqual(cf_untreated["clearance_status"], "ARRESTED_MUCUS_STASIS")
        self.assertEqual(cf_untreated["pcl_brush_compression_ratio"], 0.0)

        # Partially rescued: intermediate compromised clearance
        self.assertGreater(cf_rescued["transport_velocity_um_per_s"], 5.0)
        self.assertLess(cf_rescued["transport_velocity_um_per_s"], healthy["transport_velocity_um_per_s"])

        # PCD (immotile cilia): near zero velocity despite low mucus viscosity
        self.assertLess(pcd["transport_velocity_um_per_s"], 3.0)

    def test_viscoelastic_modulus_scaling(self) -> None:
        healthy = self.report["profiles"]["Healthy_Hydrated"]
        cf = self.report["profiles"]["CF_Untreated_Hyperconcentrated"]

        # CF complex shear modulus is > 10x higher than healthy
        self.assertGreater(cf["complex_modulus_g_star_pa"], healthy["complex_modulus_g_star_pa"] * 10.0)
        self.assertGreater(cf["viscoelastic_resistance_n_s_per_m"], healthy["viscoelastic_resistance_n_s_per_m"] * 10.0)

    def test_ah10_dynamic_falsification_negative_control(self) -> None:
        fals_result = falsify_mucociliary_clearance()
        self.assertEqual(fals_result["rule"], "AH-10")
        self.assertEqual(fals_result["verdict"], "PASS")
        self.assertTrue(fals_result["negative_control_passed"])
        self.assertEqual(fals_result["arrested_transport_velocity"], 0.0)


if __name__ == "__main__":
    unittest.main()
