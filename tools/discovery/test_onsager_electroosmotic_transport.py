"""Unit tests for Round 24 (T-0546): Onsager Reciprocal Relations & Cross-Coupled Electro-Osmotic Transport Across Airway Surface Liquid.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRev.37.405 (Reciprocal Relations in Irreversible Processes. I.)
  - DOI:10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)
  - PDB 6MSM (Human phosphorylated ATP-bound CFTR cryo-EM structure)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_onsager_electroosmotic_transport,
)


class OnsagerElectroosmoticTransportTests(unittest.TestCase):
    """Test suite verifying Round 24 Onsager reciprocal relations and ASL electro-osmotic coupling."""

    def setUp(self) -> None:
        self.report = evaluate_onsager_electroosmotic_transport()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0546")
        self.assertTrue(self.report["onsager_transport_bridge_verified"])
        self.assertTrue(self.report["epithelial_transport_stratification_verified"])
        self.assertTrue(self.report["proof_reciprocal_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_onsager_reciprocal_symmetry_across_regimes(self) -> None:
        profiles = self.report["transport_profiles"]
        self.assertIn("healthy_airway_cftr_wt", profiles)
        self.assertIn("cf_epithelium_f508del", profiles)
        self.assertIn("rescued_epithelium_trikafta", profiles)
        self.assertIn("passive_equilibrium_control", profiles)

        for reg_id, prof in profiles.items():
            self.assertTrue(
                prof["onsager_reciprocal_satisfied"],
                f"Onsager reciprocal symmetry violated for regime {reg_id}",
            )
            self.assertLessEqual(
                prof["symmetry_error"],
                1e-4,
                f"Symmetry residual exceeds 1e-4 for regime {reg_id}",
            )
            self.assertTrue(
                prof["positive_definite"],
                f"Onsager matrix is not positive-definite for regime {reg_id}",
            )
            self.assertTrue(
                prof["second_law_satisfied"],
                f"Entropy production violates Second Law for regime {reg_id}",
            )

    def test_biophysical_genotype_stratification(self) -> None:
        profiles = self.report["transport_profiles"]
        wt = profiles["healthy_airway_cftr_wt"]
        cf = profiles["cf_epithelium_f508del"]
        rescued = profiles["rescued_epithelium_trikafta"]

        # Healthy: High CFTR conductance, robust PCL hydration and coupling
        self.assertGreaterEqual(wt["hydraulic_permeability_l11"], 2.0)
        self.assertGreaterEqual(wt["electrical_conductance_l22"], 4.0)
        self.assertGreaterEqual(wt["pcl_depth_um"], 7.0)
        self.assertGreaterEqual(wt["hydration_stability_index"], 0.70)
        self.assertGreaterEqual(wt["coupling_degree_q"], 0.30)

        # CF: Collapsed conductance, severe PCL dehydration
        self.assertLess(cf["hydraulic_permeability_l11"], 1.0)
        self.assertLess(cf["electrical_conductance_l22"], 2.0)
        self.assertLess(cf["pcl_depth_um"], 3.0)
        self.assertLess(cf["hydration_stability_index"], 0.30)
        self.assertLess(cf["coupling_degree_q"], 0.25)

        # Rescued: Normalized conductance and restored PCL depth
        self.assertGreaterEqual(rescued["hydraulic_permeability_l11"], 2.0)
        self.assertGreaterEqual(rescued["electrical_conductance_l22"], 3.5)
        self.assertGreaterEqual(rescued["pcl_depth_um"], 6.5)
        self.assertGreaterEqual(rescued["hydration_stability_index"], 0.70)
        self.assertGreaterEqual(rescued["coupling_degree_q"], 0.30)

    def test_formal_proof_reciprocal_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["onsager_reciprocal_satisfied"])
            self.assertLessEqual(c_prof["symmetry_error"], 1e-4)
            self.assertGreater(c_prof["determinant"], 0.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_onsager_electroosmotic_transport(n_transport_regimes=2)


if __name__ == "__main__":
    unittest.main()
