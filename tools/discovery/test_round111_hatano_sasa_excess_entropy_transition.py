"""Unit tests for Round 111: Hatano-Sasa NESS Excess Entropy Production & Housekeeping Heat Decomposition."""

import unittest
from tools.discovery.round111_hatano_sasa_excess_entropy_transition import (
    CITATIONS,
    HYPOTHESIS,
    HatanoSasaTransitionProfile,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    integrate_exponential_fluctuation_average,
    run_falsification,
    run_falsification_checks,
    simulate_hatano_sasa_transition,
)


class TestRound111HatanoSasaExcessEntropyTransition(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round111")
        self.assertEqual(TASK_ID, "T-0743")
        self.assertIn("Hatano-Sasa", TITLE)
        self.assertIn("Housekeeping Heat", TITLE)
        self.assertIn("Hatano-Sasa", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        expected_dois = {
            "10.1103/PhysRevLett.86.3463",
            "10.1088/0305-4470/38/34/L03",
            "10.1152/physrev.00007.2018",
            "10.1085/jgp.200609558",
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
            "wt_cftr_atp_ivacaftor",
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

    def test_05_hatano_sasa_and_speck_seifert_integral_fluctuation_theorems(self) -> None:
        for p in PROFILES:
            res = simulate_hatano_sasa_transition(p)
            self.assertAlmostEqual(res.hatano_sasa_exp_average, 1.0, places=6)
            self.assertAlmostEqual(res.speck_seifert_hk_exp_average, 1.0, places=6)
            self.assertTrue(res.is_hatano_sasa_ift_valid)
            self.assertTrue(res.is_housekeeping_positivity_satisfied)

    def test_06_generalized_clausius_and_three_faces_decomposition(self) -> None:
        for p in PROFILES:
            res = simulate_hatano_sasa_transition(p)
            # Delta_S_sys >= -beta * <Q_ex>  <=>  Delta_S_sys + beta * <Q_ex> = <Y> >= 0
            self.assertGreaterEqual(res.delta_s_sys_kb, -res.excess_heat_kb_t)
            self.assertAlmostEqual(
                res.delta_s_sys_kb + res.excess_heat_kb_t,
                res.excess_entropy_y_kb,
                places=12,
            )
            self.assertAlmostEqual(
                res.housekeeping_heat_kb_t + res.excess_entropy_y_kb,
                res.total_entropy_production_kb,
                places=12,
            )

    def test_07_modulator_rescue_suppresses_excess_transition_entropy(self) -> None:
        results = {p.name: simulate_hatano_sasa_transition(p) for p in PROFILES}
        wt = results["wt_cftr_atp_ivacaftor"]
        eti = results["f508del_eti_rescued"]
        g551d = results["g551d_ivacaftor_potentiated"]
        uncorr = results["f508del_uncorrected"]
        w1282x = results["w1282x_nonsense_truncated"]

        self.assertLess(wt.excess_entropy_y_kb, eti.excess_entropy_y_kb)
        self.assertLess(eti.excess_entropy_y_kb, g551d.excess_entropy_y_kb)
        self.assertLess(g551d.excess_entropy_y_kb, uncorr.excess_entropy_y_kb)
        self.assertLess(uncorr.excess_entropy_y_kb, w1282x.excess_entropy_y_kb)
        self.assertGreater(uncorr.excess_entropy_y_kb / eti.excess_entropy_y_kb, 5.0)
        self.assertGreater(eti.adiabatic_housekeeping_fraction_pct, 94.0)

    def test_08_poisoned_controls_fail_strictly(self) -> None:
        p_bad_hs = HatanoSasaTransitionProfile(
            name="bad_hs",
            description="Negative excess entropy violation",
            q_hk_rate_kb_t_s=15.0,
            transition_tau_s=0.50,
            sigma_y_kb=1.0,
            delta_s_sys_kb=-1.0,
            poison_mode="broken_hatano_sasa_ift",
        )
        res_hs = simulate_hatano_sasa_transition(p_bad_hs)
        self.assertFalse(res_hs.is_generalized_clausius_satisfied)
        self.assertFalse(res_hs.is_hatano_sasa_ift_valid)
        self.assertLess(res_hs.excess_entropy_y_kb, 0.0)

        p_bad_hk = HatanoSasaTransitionProfile(
            name="bad_hk",
            description="Negative housekeeping heat violation",
            q_hk_rate_kb_t_s=-5.0,
            transition_tau_s=0.50,
            sigma_y_kb=1.0,
            delta_s_sys_kb=-1.0,
            poison_mode="negative_housekeeping_heat",
        )
        res_hk = simulate_hatano_sasa_transition(p_bad_hk)
        self.assertFalse(res_hk.is_housekeeping_positivity_satisfied)
        self.assertLess(res_hk.housekeeping_heat_kb_t, 0.0)
        self.assertGreater(integrate_exponential_fluctuation_average(-1.0, 1.0), 2.0)


if __name__ == "__main__":
    unittest.main()
