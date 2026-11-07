"""Unit tests for Round 110: Gallavotti-Cohen Stochastic Action & Lebowitz-Spohn Fluctuation Symmetry."""

import math
import unittest
from tools.discovery.round110_gallavotti_cohen_stochastic_action import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    StochasticActionProfile,
    TASK_ID,
    TITLE,
    evaluate_large_deviation_rate_function,
    evaluate_lebowitz_spohn_scgf,
    run_falsification,
    run_falsification_checks,
    simulate_stochastic_action,
)


class TestRound110GallavottiCohenStochasticAction(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round110")
        self.assertEqual(TASK_ID, "T-0741")
        self.assertIn("Gallavotti-Cohen", TITLE)
        self.assertIn("Lebowitz-Spohn", TITLE)
        self.assertIn("Lebowitz-Spohn", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        expected_dois = {
            "10.1023/A:1004589714161",
            "10.1088/0034-4885/75/12/126001",
            "10.1085/jgp.200609558",
            "10.1038/nphys3230",
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
        checks = report["checks"]
        for check_name, check_val in checks.items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

    def test_05_lebowitz_spohn_scgf_exact_involution(self) -> None:
        for p in PROFILES:
            affinity = math.log(p.j_plus_s_inv / p.j_minus_s_inv)
            for k in (-1.5, -0.7, 0.0, 0.4, 1.1):
                lhs = evaluate_lebowitz_spohn_scgf(k, p.j_plus_s_inv, p.j_minus_s_inv)
                rhs = evaluate_lebowitz_spohn_scgf(-affinity - k, p.j_plus_s_inv, p.j_minus_s_inv)
                self.assertAlmostEqual(lhs, rhs, places=10)

    def test_06_large_deviation_rate_function_fluctuation_relation(self) -> None:
        for p in PROFILES:
            affinity = math.log(p.j_plus_s_inv / p.j_minus_s_inv)
            mean_j = p.j_plus_s_inv - p.j_minus_s_inv
            # Rate function vanishes at <J>
            i_mean = evaluate_large_deviation_rate_function(mean_j, p.j_plus_s_inv, p.j_minus_s_inv)
            self.assertAlmostEqual(i_mean, 0.0, places=12)
            # I(-j) - I(j) = A * j
            for j_val in (0.25 * mean_j, 0.8 * mean_j, 1.6 * mean_j):
                i_pos = evaluate_large_deviation_rate_function(j_val, p.j_plus_s_inv, p.j_minus_s_inv)
                i_neg = evaluate_large_deviation_rate_function(-j_val, p.j_plus_s_inv, p.j_minus_s_inv)
                self.assertAlmostEqual(i_neg - i_pos, affinity * j_val, places=10)

    def test_07_zero_current_stalling_barrier_identity(self) -> None:
        # At j = 0, I(0) = (sqrt(j_+) - sqrt(j_-))^2
        for p in PROFILES:
            res = simulate_stochastic_action(p)
            expected_i0 = (math.sqrt(p.j_plus_s_inv) - math.sqrt(p.j_minus_s_inv)) ** 2
            self.assertAlmostEqual(res.rate_function_at_zero_current, expected_i0, places=10)

    def test_08_modulator_rescue_and_potentiation_hierarchy(self) -> None:
        results = {p.name: simulate_stochastic_action(p) for p in PROFILES}
        wt = results["wt_cftr_atp_driven"]
        eti = results["f508del_eti_rescued"]
        g551d = results["g551d_ivacaftor_potentiated"]
        uncorr = results["f508del_uncorrected"]
        w1282x = results["w1282x_nonsense_truncated"]

        self.assertGreater(wt.mean_current_s_inv, eti.mean_current_s_inv)
        self.assertGreater(eti.mean_current_s_inv, g551d.mean_current_s_inv)
        self.assertGreater(g551d.mean_current_s_inv, uncorr.mean_current_s_inv)
        self.assertGreater(uncorr.mean_current_s_inv, w1282x.mean_current_s_inv)
        self.assertGreater(eti.mean_current_s_inv / uncorr.mean_current_s_inv, 10.0)

    def test_09_poisoned_broken_gc_symmetry_rejected(self) -> None:
        p_bad_gc = StochasticActionProfile(
            name="bad_gc",
            description="Broken GC symmetry",
            j_plus_s_inv=40.0,
            j_minus_s_inv=2.0,
            delta_mu_atp_kb_t=2.0,
            electrochemical_driving_mv=-26.7,
            poison_mode="broken_gc_symmetry",
        )
        res = simulate_stochastic_action(p_bad_gc)
        self.assertFalse(res.is_gc_symmetry_valid)
        self.assertGreater(res.gc_scgf_symmetry_max_error, 1.0)
        self.assertGreater(res.gc_rate_function_symmetry_max_error, 0.1)

    def test_10_poisoned_negative_affinity_second_law_rejected(self) -> None:
        p_bad_aff = StochasticActionProfile(
            name="bad_affinity",
            description="Negative affinity violation",
            j_plus_s_inv=0.5,
            j_minus_s_inv=5.0,
            delta_mu_atp_kb_t=-1.5,
            electrochemical_driving_mv=20.0,
            poison_mode="negative_affinity",
        )
        res = simulate_stochastic_action(p_bad_aff)
        self.assertFalse(res.is_second_law_satisfied)
        self.assertLess(res.affinity_kb_t, 0.0)
        self.assertLess(res.entropy_production_rate_kb_s, 0.0)


if __name__ == "__main__":
    unittest.main()
