"""Unit tests for Round 116: Cross-Field World Economic Forum Global Resilience Bridge — Thermostable LNP Pareto & Equity."""

import unittest
from tools.discovery.round116_wef_global_health_coldchain_pareto_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    TOTAL_GLOBAL_BATCH_DOSES,
    WefColdChainProfile,
    compute_arrhenius_degradation_rate,
    compute_gini_coefficient,
    run_falsification,
    run_falsification_checks,
    run_round116_verification,
    simulate_wef_coldchain_allocation,
)


class TestRound116WefGlobalHealthColdchainParetoBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round116")
        self.assertEqual(TASK_ID, "T-0751")
        self.assertEqual(PRIMARY_DOI, "10.1021/acs.nanolett.6b04242")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1016/j.bioactmat.2020.03.001", "10.1056/NEJMms2025396"],
        )
        self.assertIn("World Economic Forum", TITLE)
        self.assertIn("Pareto-Optimal", TITLE)
        self.assertIn("Arrhenius", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 3)
        expected_dois = {
            "10.1021/acs.nanolett.6b04242",
            "10.1016/j.bioactmat.2020.03.001",
            "10.1056/NEJMms2025396",
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

    def test_03_five_wef_coldchain_profiles_inventory(self) -> None:
        expected_names = [
            "lyophilized_ambient_equitable_pareto",
            "regional_hub_dual_chain_eti_lnp",
            "standard_minus20c_regional_network",
            "ultracold_minus80c_centralized_bottleneck",
            "unbuffered_liquid_lnp_spoilage",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round116_verification()
        self.assertEqual(ver["round"], 116)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1021/acs.nanolett.6b04242")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_arrhenius_kinetics_mass_conservation_and_gini_bounds(self) -> None:
        for p in PROFILES:
            res = simulate_wef_coldchain_allocation(p)
            self.assertGreater(res.arrhenius_degradation_rate_day_inv, 0.0)
            self.assertGreater(res.intact_bioactive_fraction, 0.0)
            self.assertLessEqual(res.intact_bioactive_fraction, 1.0)
            self.assertAlmostEqual(res.total_allocated_doses, TOTAL_GLOBAL_BATCH_DOSES, places=9)
            self.assertLess(res.mass_conservation_error_doses, 1e-9)
            self.assertGreaterEqual(res.gini_inequality_index, 0.0)
            self.assertLessEqual(res.gini_inequality_index, 1.0)
            self.assertAlmostEqual(
                res.gini_inequality_index + res.wef_health_equity_index,
                1.0,
                places=12,
            )
            self.assertTrue(res.is_arrhenius_stability_valid)
            self.assertTrue(res.is_gini_equity_valid)
            self.assertTrue(res.is_mass_conservation_satisfied)

    def test_06_lyophilized_pareto_dominance_and_global_equity_ordering(self) -> None:
        results = {p.name: simulate_wef_coldchain_allocation(p) for p in PROFILES}
        lyo = results["lyophilized_ambient_equitable_pareto"]
        hub = results["regional_hub_dual_chain_eti_lnp"]
        minus20 = results["standard_minus20c_regional_network"]
        minus80 = results["ultracold_minus80c_centralized_bottleneck"]
        spoiled = results["unbuffered_liquid_lnp_spoilage"]

        self.assertGreater(lyo.wef_social_welfare_score, hub.wef_social_welfare_score)
        self.assertGreater(hub.wef_social_welfare_score, minus20.wef_social_welfare_score)
        self.assertGreater(minus20.wef_social_welfare_score, minus80.wef_social_welfare_score)
        self.assertGreater(minus80.wef_social_welfare_score, spoiled.wef_social_welfare_score)

        self.assertGreater(
            lyo.global_potent_coverage_fraction / minus80.global_potent_coverage_fraction,
            3.5,
        )
        self.assertGreater(lyo.wef_health_equity_index, 0.92)
        self.assertLess(minus80.wef_health_equity_index, 0.45)

    def test_07_three_poisoned_negative_controls_fail_strictly(self) -> None:
        p_bio = WefColdChainProfile(
            name="bad_bio",
            description="Super-unitary bioactive survival violation",
            activation_energy_ea_kj_mol=90.0,
            arrhenius_pre_exp_day_inv=8.0e13,
            transit_temp_kelvin=293.15,
            transit_days=30.0,
            regional_dose_shares=(0.25, 0.25, 0.25, 0.25),
            regional_last_mile_reach=(0.9, 0.9, 0.9, 0.9),
            cost_carbon_efficiency=0.8,
            poison_mode="super_unitary_bioactive_fraction",
        )
        res_bio = simulate_wef_coldchain_allocation(p_bio)
        self.assertFalse(res_bio.is_arrhenius_stability_valid)
        self.assertGreater(res_bio.intact_bioactive_fraction, 1.0)

        p_gini = WefColdChainProfile(
            name="bad_gini",
            description="Invalid Gini bound violation",
            activation_energy_ea_kj_mol=90.0,
            arrhenius_pre_exp_day_inv=8.0e13,
            transit_temp_kelvin=293.15,
            transit_days=30.0,
            regional_dose_shares=(0.25, 0.25, 0.25, 0.25),
            regional_last_mile_reach=(0.9, 0.9, 0.9, 0.9),
            cost_carbon_efficiency=0.8,
            poison_mode="invalid_gini_equity_bound",
        )
        res_gini = simulate_wef_coldchain_allocation(p_gini)
        self.assertFalse(res_gini.is_gini_equity_valid)
        self.assertGreater(res_gini.wef_health_equity_index, 1.0)

        p_mass = WefColdChainProfile(
            name="bad_mass",
            description="Broken mass conservation violation",
            activation_energy_ea_kj_mol=90.0,
            arrhenius_pre_exp_day_inv=8.0e13,
            transit_temp_kelvin=293.15,
            transit_days=30.0,
            regional_dose_shares=(0.25, 0.25, 0.25, 0.25),
            regional_last_mile_reach=(0.9, 0.9, 0.9, 0.9),
            cost_carbon_efficiency=0.8,
            poison_mode="broken_mass_conservation",
        )
        res_mass = simulate_wef_coldchain_allocation(p_mass)
        self.assertFalse(res_mass.is_mass_conservation_satisfied)
        self.assertGreater(res_mass.mass_conservation_error_doses, 1000.0)
        self.assertAlmostEqual(compute_gini_coefficient((100.0, 100.0, 100.0, 100.0)), 0.0, places=12)
        self.assertGreater(compute_arrhenius_degradation_rate(80.0, 1.0e12, 300.0), 0.0)


if __name__ == "__main__":
    unittest.main()
