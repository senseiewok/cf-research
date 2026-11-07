"""Unit tests for Round 113: Esposito-Van den Broeck Three-Way Entropy Production Decomposition."""

import unittest
from tools.discovery.round113_esposito_van_den_broeck_three_faces import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    EspositoThreeFacesProfile,
    integrate_ift_fluctuation_average,
    run_falsification,
    run_falsification_checks,
    simulate_esposito_three_faces,
)


class TestRound113EspositoVanDenBroeckThreeFaces(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round113")
        self.assertEqual(TASK_ID, "T-0746")
        self.assertIn("Esposito-Van den Broeck", TITLE)
        self.assertIn("Housekeeping Heat", TITLE)
        self.assertIn("Esposito-Van den Broeck", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        expected_dois = {
            "10.1103/PhysRevLett.104.090601",
            "10.1103/PhysRevE.82.011143",
            "10.1085/jgp.200609558",
            "10.1152/physrev.00007.2018",
            "10.1126/science.ade2216",
        }
        self.assertEqual({c["doi"] for c in CITATIONS}, expected_dois)
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

    def test_03_five_cftr_profiles_inventory(self) -> None:
        expected_names = [
            "wt_cftr_atp_driven",
            "f508del_eti_rescued",
            "g551d_ivacaftor_potentiated",
            "f508del_uncorrected",
            "w1282x_nonsense_truncated",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_pass(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

    def test_05_esposito_identity_and_second_law(self) -> None:
        for p in PROFILES:
            res = simulate_esposito_three_faces(p)
            self.assertLess(res.esposito_identity_residual_kb_s, 1e-12)
            self.assertTrue(res.is_esposito_identity_valid)
            self.assertTrue(res.is_second_law_satisfied)
            self.assertGreaterEqual(res.housekeeping_entropy_rate_kb_s, 0.0)
            self.assertGreaterEqual(res.excess_entropy_rate_kb_s, 0.0)
            self.assertGreaterEqual(res.total_entropy_rate_kb_s, 0.0)
            self.assertAlmostEqual(
                res.housekeeping_heat_power_zw + res.excess_heat_power_zw,
                res.total_dissipation_power_zw,
                places=3,
            )

    def test_06_three_independent_integral_fluctuation_theorems(self) -> None:
        for p in PROFILES:
            res = simulate_esposito_three_faces(p)
            self.assertTrue(res.is_all_three_ift_valid)
            self.assertAlmostEqual(res.total_ift_exp_average, 1.0, places=4)
            self.assertAlmostEqual(res.excess_ift_exp_average, 1.0, places=4)
            self.assertAlmostEqual(res.housekeeping_ift_exp_average, 1.0, places=4)

    def test_07_eti_rescue_restores_ness_housekeeping_efficiency(self) -> None:
        results = {p.name: simulate_esposito_three_faces(p) for p in PROFILES}
        wt = results["wt_cftr_atp_driven"]
        eti = results["f508del_eti_rescued"]
        g551d = results["g551d_ivacaftor_potentiated"]
        uncorr = results["f508del_uncorrected"]
        w1282x = results["w1282x_nonsense_truncated"]

        # WT and rescued profiles channel most energy into housekeeping flux
        self.assertGreater(wt.housekeeping_efficiency_pct, 88.0)
        self.assertGreater(eti.housekeeping_efficiency_pct, 80.0)
        self.assertGreater(g551d.housekeeping_efficiency_pct, 70.0)
        self.assertLess(uncorr.housekeeping_efficiency_pct, 40.0)
        self.assertLess(w1282x.housekeeping_efficiency_pct, 20.0)

        # Futile excess dissipation is suppressed in rescued states
        self.assertLess(eti.excess_entropy_rate_kb_s, uncorr.excess_entropy_rate_kb_s)
        self.assertGreater(
            uncorr.excess_entropy_rate_kb_s / eti.excess_entropy_rate_kb_s, 1.5
        )

        # Housekeeping power ranking reflects catalytic competence
        self.assertGreater(wt.housekeeping_heat_power_zw, eti.housekeeping_heat_power_zw)
        self.assertGreater(eti.housekeeping_heat_power_zw, g551d.housekeeping_heat_power_zw)
        self.assertGreater(g551d.housekeeping_heat_power_zw, uncorr.housekeeping_heat_power_zw)
        self.assertGreater(uncorr.housekeeping_heat_power_zw, w1282x.housekeeping_heat_power_zw)


if __name__ == "__main__":
    unittest.main()
