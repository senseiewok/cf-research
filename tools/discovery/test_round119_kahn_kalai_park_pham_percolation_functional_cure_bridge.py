"""Unit tests for Round 119: Park-Pham / Kahn-Kalai Hypergraph Expectation Threshold & CF Epithelial Percolation Cure Boundary."""

import unittest
from tools.discovery.round119_kahn_kalai_park_pham_percolation_functional_cure_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    KahnKalaiPercolationProfile,
    run_falsification,
    run_falsification_checks,
    run_round119_verification,
    simulate_kahn_kalai_percolation,
)


class TestRound119KahnKalaiParkPhamPercolationFunctionalCureBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round119")
        self.assertEqual(TASK_ID, "T-0755")
        self.assertEqual(PRIMARY_DOI, "10.4007/annals.2024.199.1.4")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1038/ng0992-21", "10.1164/rccm.202308-1363OC"],
        )
        self.assertIn("Park-Pham", TITLE)
        self.assertIn("11.8%", TITLE)
        self.assertIn("Kahn-Kalai", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 3)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertTrue(ah10["passed"])
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], ROUND_ID)
        self.assertEqual(ah10["task_id"], TASK_ID)

    def test_02_parameter_provenance_bounds(self) -> None:
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn("unit", p_meta, msg=f"Missing unit for {p_name}")
            self.assertIn("kind", p_meta, msg=f"Missing kind for {p_name}")
            self.assertIn("source", p_meta, msg=f"Missing source for {p_name}")
            self.assertIn("bounds", p_meta, msg=f"Missing bounds for {p_name}")
            low, high = p_meta["bounds"]
            self.assertLess(low, high)
            if "value" in p_meta:
                val = p_meta["value"]
                self.assertGreaterEqual(val, low)
                self.assertLessEqual(val, high)

    def test_03_five_percolation_profiles_inventory(self) -> None:
        expected_names = [
            "uniform_peg_lnp_basal_ciliated_cure",
            "homogeneous_muco_inert_mrna_above_threshold",
            "critical_percolation_transition_regimen",
            "patchy_proximal_lnp_sub_percolation_trap",
            "unshielded_vector_focal_failure",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round119_verification()
        self.assertEqual(ver["round"], 119)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.4007/annals.2024.199.1.4")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_park_pham_expectation_threshold_and_16pct_uniform_vs_patchy_contrast(self) -> None:
        results = {p.name: simulate_kahn_kalai_percolation(p) for p in PROFILES}
        homog = results["homogeneous_muco_inert_mrna_above_threshold"]
        patchy = results["patchy_proximal_lnp_sub_percolation_trap"]

        self.assertAlmostEqual(homog.critical_percolation_pc, 0.1180, places=4)
        self.assertAlmostEqual(homog.transfected_cell_fraction, patchy.transfected_cell_fraction, places=12)
        self.assertGreater(homog.hypergraph_percolation_probability, 0.95)
        self.assertLess(patchy.hypergraph_percolation_probability, 0.06)

    def test_06_monotonic_ordering_and_negative_controls(self) -> None:
        results = [simulate_kahn_kalai_percolation(p) for p in PROFILES]
        for i in range(len(results) - 1):
            self.assertGreater(
                results[i].effective_percolating_fraction,
                results[i + 1].effective_percolating_fraction,
            )
            self.assertGreater(
                results[i].hypergraph_percolation_probability,
                results[i + 1].hypergraph_percolation_probability,
            )
            self.assertGreater(
                results[i].restored_isc_pct_wt,
                results[i + 1].restored_isc_pct_wt,
            )


if __name__ == "__main__":
    unittest.main()
