"""Unit tests for Round 109: Crooks Non-Equilibrium Work Fluctuation Theorem & Jarzynski Free-Energy Recovery."""

import unittest
from tools.discovery.round109_crooks_jarzynski_nonequilibrium_folding import (
    CITATIONS,
    HYPOTHESIS,
    FluctuationStateProfile,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    run_falsification,
    run_falsification_checks,
    simulate_nonequilibrium_work,
)


class TestRound109CrooksJarzynskiNonequilibriumFolding(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round109")
        self.assertEqual(TASK_ID, "T-0740")
        self.assertIn("Crooks", TITLE)
        self.assertIn("Jarzynski", TITLE)
        self.assertIn("NBD1-NBD2", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertTrue(ah10["passed"])
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], ROUND_ID)
        self.assertEqual(ah10["task_id"], TASK_ID)

    def test_02_parameter_provenance_bounds(self) -> None:
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn("unit", p_meta)
            self.assertIn("kind", p_meta)
            self.assertIn("source", p_meta)
            self.assertIn("bounds", p_meta)
            low, high = p_meta["bounds"]
            self.assertLess(low, high)
            if "value" in p_meta:
                val = p_meta["value"]
                self.assertGreaterEqual(val, low)
                self.assertLessEqual(val, high)

    def test_03_all_checks_pass(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        checks = report["checks"]
        self.assertTrue(checks["all_valid_profiles_satisfy_second_law"])
        self.assertTrue(checks["all_valid_profiles_satisfy_crooks_symmetry"])
        self.assertTrue(checks["jarzynski_recovers_free_energy_within_1pct"])
        self.assertTrue(checks["crooks_crossing_matches_equilibrium_free_energy"])
        self.assertTrue(checks["wt_potentiator_suppresses_dissipation_below_1_kb_t"])
        self.assertTrue(checks["trikafta_restores_dissipation_over_2x_relative_to_f508del"])
        self.assertTrue(checks["f508del_dissipation_exceeds_8_kb_t"])
        self.assertTrue(checks["poisoned_perpetual_motion_violates_second_law"])
        self.assertTrue(checks["poisoned_broken_crooks_violates_microreversibility"])
        self.assertTrue(checks["parametric_sensitivity_25pct_preserves_ordering"])

    def test_04_individual_simulation_accuracy(self) -> None:
        wt = simulate_nonequilibrium_work(PROFILES[0])
        self.assertEqual(wt.profile_name, "wildtype_atp_potentiator")
        self.assertTrue(wt.is_second_law_satisfied)
        self.assertTrue(wt.is_crooks_symmetry_valid)
        self.assertLess(wt.dissipated_work_kb_t, 1.0)
        self.assertLess(wt.jarzynski_error_pct, 0.5)
        self.assertLess(wt.crooks_symmetry_error_pct, 0.05)
        self.assertAlmostEqual(wt.fdt_work_variance_ratio, 1.0, places=6)

    def test_05_poisoned_controls_fail_strictly(self) -> None:
        p_bad = FluctuationStateProfile(
            name="bad_control",
            description="Perpetual motion control",
            delta_f_eq_kcal_mol=-10.0,
            sigma_w_kb_t=1.0,
            driving_protocol_speed_ratio=1.0,
            is_poisoned=True,
        )
        res_bad = simulate_nonequilibrium_work(p_bad)
        self.assertFalse(res_bad.is_second_law_satisfied)
        self.assertFalse(res_bad.is_crooks_symmetry_valid)
        self.assertLess(res_bad.dissipated_work_kb_t, 0.0)


if __name__ == "__main__":
    unittest.main()
