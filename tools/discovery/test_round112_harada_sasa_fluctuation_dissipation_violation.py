"""Unit tests for Round 112: Harada-Sasa Spectral Fluctuation-Response Violation Integral & Active ATP Power."""

import unittest
from tools.discovery.round112_harada_sasa_fluctuation_dissipation_violation import (
    CITATIONS,
    HYPOTHESIS,
    HaradaSasaGatingProfile,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    evaluate_spectral_fdt_violation_density,
    integrate_harada_sasa_spectrum,
    run_falsification,
    run_falsification_checks,
    simulate_harada_sasa_gating,
)


class TestRound112HaradaSasaFluctuationDissipationViolation(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round112")
        self.assertEqual(TASK_ID, "T-0744")
        self.assertIn("Harada-Sasa", TITLE)
        self.assertIn("Fluctuation-Response", TITLE)
        self.assertIn("Harada-Sasa", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 5)
        expected_dois = {
            "10.1103/PhysRevLett.95.130602",
            "10.1103/PhysRevLett.104.198103",
            "10.1038/nature03313",
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

    def test_05_harada_sasa_spectral_quadrature_and_power_balance(self) -> None:
        for p in PROFILES:
            res = simulate_harada_sasa_gating(p)
            self.assertLess(res.spectral_quadrature_error, 1e-12)
            self.assertLess(res.harada_sasa_power_balance_error, 1e-12)
            self.assertTrue(res.is_harada_sasa_equality_valid)
            self.assertTrue(res.is_second_law_satisfied)
            self.assertAlmostEqual(
                res.directed_drift_power_kbt_s + res.spectral_dissipation_power_kbt_s,
                res.total_harada_sasa_dissipation_kbt_s,
                places=12,
            )

    def test_06_modulator_rescue_restores_directed_drift_and_stokes_efficiency(self) -> None:
        results = {p.name: simulate_harada_sasa_gating(p) for p in PROFILES}
        wt = results["wt_cftr_atp_driven"]
        eti = results["f508del_eti_rescued"]
        g551d = results["g551d_ivacaftor_potentiated"]
        uncorr = results["f508del_uncorrected"]
        w1282x = results["w1282x_nonsense_truncated"]

        self.assertGreater(wt.directed_drift_power_kbt_s, eti.directed_drift_power_kbt_s)
        self.assertGreater(eti.directed_drift_power_kbt_s, g551d.directed_drift_power_kbt_s)
        self.assertGreater(g551d.directed_drift_power_kbt_s, uncorr.directed_drift_power_kbt_s)
        self.assertGreater(uncorr.directed_drift_power_kbt_s, w1282x.directed_drift_power_kbt_s)
        self.assertGreater(eti.directed_drift_power_kbt_s / uncorr.directed_drift_power_kbt_s, 15.0)
        self.assertGreater(eti.stokes_mechanochemical_efficiency_pct, 50.0)

    def test_07_poisoned_controls_fail_strictly(self) -> None:
        p_neg = HaradaSasaGatingProfile(
            name="bad_neg",
            description="Negative FDT violation",
            gamma_kbt_s_nm2=2.0,
            v_ss_nm_s=0.50,
            delta_act_kb_t=-3.20,
            omega_c_rad_s=10.0,
            poison_mode="negative_fdt_violation",
        )
        res_neg = simulate_harada_sasa_gating(p_neg)
        self.assertFalse(res_neg.is_second_law_satisfied)
        self.assertLess(res_neg.spectral_fdt_integral_nm2_s2, 0.0)

        p_broken = HaradaSasaGatingProfile(
            name="bad_balance",
            description="Broken Harada-Sasa power balance",
            gamma_kbt_s_nm2=2.0,
            v_ss_nm_s=3.00,
            delta_act_kb_t=1.50,
            omega_c_rad_s=12.0,
            poison_mode="broken_power_balance",
        )
        res_broken = simulate_harada_sasa_gating(p_broken)
        self.assertFalse(res_broken.is_harada_sasa_equality_valid)
        self.assertGreater(res_broken.harada_sasa_power_balance_error, 10.0)
        self.assertGreater(evaluate_spectral_fdt_violation_density(0.0, 2.0, 1.8, 15.0), 1.0)
        self.assertAlmostEqual(integrate_harada_sasa_spectrum(2.0, 1.8, 15.0), 13.5, places=12)


if __name__ == "__main__":
    unittest.main()
