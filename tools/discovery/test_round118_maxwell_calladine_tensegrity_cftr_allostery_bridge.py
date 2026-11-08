"""Unit tests for Round 118: Maxwell-Calladine Topological Mechanics & Tensegrity Nullspace Rigidity in CFTR Allostery."""

import unittest
from tools.discovery.round118_maxwell_calladine_tensegrity_cftr_allostery_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    MaxwellCalladineProfile,
    run_falsification,
    run_falsification_checks,
    run_round118_verification,
    simulate_maxwell_calladine_tensegrity,
)


class TestRound118MaxwellCalladineTensegrityCftrAllosteryBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round118")
        self.assertEqual(TASK_ID, "T-0754")
        self.assertEqual(PRIMARY_DOI, "10.1016/0020-7683(86)90014-4")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1038/nphys2835", "10.1126/science.ade2216"],
        )
        self.assertIn("Maxwell-Calladine", TITLE)
        self.assertIn("Tensegrity", TITLE)
        self.assertIn("N_floppy - N_self_stress", HYPOTHESIS)
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

    def test_03_five_tensegrity_profiles_inventory(self) -> None:
        expected_names = [
            "wt_cftr_atp_bound_tensegrity_locked",
            "f508del_eti_triple_strut_rescued",
            "f508del_lumacaftor_iva_dual_strut",
            "f508del_uncorrected_icl4_unpinned",
            "w1282x_truncated_nbd2_hypostatic",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round118_verification()
        self.assertEqual(ver["round"], 118)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1016/0020-7683(86)90014-4")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_maxwell_calladine_topological_index_identity(self) -> None:
        for p in PROFILES:
            res = simulate_maxwell_calladine_tensegrity(p)
            self.assertEqual(
                res.floppy_modes_count - res.self_stress_modes_count,
                3 * res.vertices_v - res.struts_e - 6,
            )
            self.assertGreaterEqual(res.floppy_modes_count, 0)
            self.assertGreaterEqual(res.self_stress_modes_count, 0)
            self.assertTrue(res.is_maxwell_calladine_satisfied)
            self.assertTrue(res.is_nullspace_non_negative)
            self.assertTrue(res.is_transmission_bound_valid)

    def test_06_monotonic_ordering_and_negative_controls(self) -> None:
        results = {p.name: simulate_maxwell_calladine_tensegrity(p) for p in PROFILES}
        wt = results["wt_cftr_atp_bound_tensegrity_locked"]
        eti = results["f508del_eti_triple_strut_rescued"]
        dual = results["f508del_lumacaftor_iva_dual_strut"]
        uncorr = results["f508del_uncorrected_icl4_unpinned"]
        w1282x = results["w1282x_truncated_nbd2_hypostatic"]

        self.assertLess(wt.floppy_modes_count, eti.floppy_modes_count)
        self.assertLess(eti.floppy_modes_count, dual.floppy_modes_count)
        self.assertLess(dual.floppy_modes_count, uncorr.floppy_modes_count)
        self.assertLess(uncorr.floppy_modes_count, w1282x.floppy_modes_count)

        self.assertGreater(wt.allosteric_transmission_efficiency, eti.allosteric_transmission_efficiency)
        self.assertGreater(eti.allosteric_transmission_efficiency, dual.allosteric_transmission_efficiency)
        self.assertGreater(dual.allosteric_transmission_efficiency, uncorr.allosteric_transmission_efficiency)
        self.assertGreater(uncorr.allosteric_transmission_efficiency, w1282x.allosteric_transmission_efficiency)


if __name__ == "__main__":
    unittest.main()
