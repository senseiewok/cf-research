"""Unit tests for Round 120: Kuramoto Ciliary Metachronal Phase-Locking & Ames Control Barrier Function Safe Dosing."""

import unittest
from tools.discovery.round120_kuramoto_cbf_ciliary_metachronal_safety_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    KuramotoCbfProfile,
    run_falsification,
    run_falsification_checks,
    run_round120_verification,
    simulate_kuramoto_cbf_bridge,
)


class TestRound120KuramotoCbfCiliaryMetachronalSafetyBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round120")
        self.assertEqual(TASK_ID, "T-0756")
        self.assertEqual(PRIMARY_DOI, "10.1038/s41567-019-0770-2")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1109/TAC.2016.2638961", "10.1056/NEJMoa1908639"],
        )
        self.assertIn("Kuramoto", TITLE)
        self.assertIn("Control Barrier Function", TITLE)
        self.assertIn("R_Kuramoto", HYPOTHESIS)
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

    def test_03_five_kuramoto_cbf_profiles_inventory(self) -> None:
        expected_names = [
            "wt_synchronized_metachronal_wave",
            "class1_cbf_qp_safe_repeat_lnp_metachronal",
            "f508del_eti_rescued_ciliary_locking",
            "partially_coupled_suboptimal_redosing",
            "dehydrated_cf_ciliary_quiescence_unlocked",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round120_verification()
        self.assertEqual(ver["round"], 120)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1038/s41567-019-0770-2")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_kuramoto_phase_transition_and_cbf_qp_safety_invariance(self) -> None:
        results = {p.name: simulate_kuramoto_cbf_bridge(p) for p in PROFILES}
        cbf_qp = results["class1_cbf_qp_safe_repeat_lnp_metachronal"]
        unlocked = results["dehydrated_cf_ciliary_quiescence_unlocked"]

        self.assertLess(cbf_qp.cbf_qp_optimal_dose_u, cbf_qp.nominal_dose_u)
        self.assertAlmostEqual(cbf_qp.cbf_nagumo_derivative_lhs, 0.0, places=12)
        self.assertEqual(unlocked.kuramoto_order_parameter, 0.0)
        self.assertEqual(unlocked.mucociliary_clearance_velocity_um_s, 0.0)

    def test_06_monotonic_ordering_and_negative_controls(self) -> None:
        results = [simulate_kuramoto_cbf_bridge(p) for p in PROFILES]
        for i in range(len(results) - 1):
            self.assertGreater(
                results[i].kuramoto_order_parameter,
                results[i + 1].kuramoto_order_parameter,
            )
            self.assertGreater(
                results[i].mucociliary_clearance_velocity_um_s,
                results[i + 1].mucociliary_clearance_velocity_um_s,
            )
            self.assertGreater(
                results[i].cbf_safety_barrier_margin,
                results[i + 1].cbf_safety_barrier_margin,
            )


if __name__ == "__main__":
    unittest.main()
