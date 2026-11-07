"""Unit tests for Round 114: Sagawa-Ueda Generalized Second Law with Information Feedback & Allosteric Mutual Information."""

import unittest
from tools.discovery.round114_sagawa_ueda_information_thermodynamics_feedback import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    SagawaUedaAllostericProfile,
    integrate_generalized_jarzynski_average,
    run_falsification,
    run_falsification_checks,
    run_round114_verification,
    simulate_sagawa_ueda_feedback,
)


class TestRound114SagawaUedaInformationThermodynamicsFeedback(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round114")
        self.assertEqual(TASK_ID, "T-0749")
        self.assertEqual(PRIMARY_DOI, "10.1103/PhysRevLett.104.090602")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1038/nphys1821", "10.1016/j.cell.2017.02.024"],
        )
        self.assertIn("Sagawa-Ueda", TITLE)
        self.assertIn("Information Feedback", TITLE)
        self.assertIn("Sagawa-Ueda", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 3)
        expected_dois = {
            "10.1103/PhysRevLett.104.090602",
            "10.1038/nphys1821",
            "10.1016/j.cell.2017.02.024",
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

    def test_03_five_cftr_genotypes_inventory(self) -> None:
        expected_names = [
            "wt_cftr_atp_sensing",
            "f508del_eti_rescued",
            "g551d_ivacaftor_potentiated",
            "f508del_uncorrected",
            "w1282x_nonsense_truncated",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round114_verification()
        self.assertEqual(ver["round"], 114)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1103/PhysRevLett.104.090602")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertLess(ver["summary_metrics"]["max_generalized_jarzynski_error"], 1e-6)
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_sagawa_ueda_generalized_jarzynski_and_second_law_bounds(self) -> None:
        for p in PROFILES:
            res = simulate_sagawa_ueda_feedback(p)
            self.assertGreater(res.mutual_information_nats, 0.0)
            self.assertGreaterEqual(res.generalized_second_law_slack_kb_t, 0.0)
            self.assertLessEqual(
                res.extractable_allosteric_work_kb_t,
                res.information_energy_bound_kb_t + 1e-12,
            )
            self.assertGreater(res.sagawa_ueda_efficiency, 0.0)
            self.assertLessEqual(res.sagawa_ueda_efficiency, 1.0)
            self.assertAlmostEqual(res.generalized_jarzynski_exp_average, 1.0, places=6)
            self.assertAlmostEqual(
                res.mean_dissipated_work_kb_t + res.mutual_information_nats,
                res.generalized_second_law_slack_kb_t,
                places=12,
            )
            self.assertTrue(res.is_generalized_second_law_satisfied)
            self.assertTrue(res.is_efficiency_within_unitary_bound)
            self.assertTrue(res.is_generalized_jarzynski_valid)

    def test_06_modulator_rescue_restores_nbd_tmd_mutual_information_and_work(self) -> None:
        results = {p.name: simulate_sagawa_ueda_feedback(p) for p in PROFILES}
        wt = results["wt_cftr_atp_sensing"]
        eti = results["f508del_eti_rescued"]
        g551d = results["g551d_ivacaftor_potentiated"]
        uncorr = results["f508del_uncorrected"]
        w1282x = results["w1282x_nonsense_truncated"]

        self.assertGreater(wt.mutual_information_nats, eti.mutual_information_nats)
        self.assertGreater(eti.mutual_information_nats, g551d.mutual_information_nats)
        self.assertGreater(g551d.mutual_information_nats, uncorr.mutual_information_nats)
        self.assertGreater(uncorr.mutual_information_nats, w1282x.mutual_information_nats)

        self.assertGreater(
            wt.extractable_allosteric_work_kb_t, eti.extractable_allosteric_work_kb_t
        )
        self.assertGreater(
            eti.extractable_allosteric_work_kb_t, g551d.extractable_allosteric_work_kb_t
        )
        self.assertGreater(
            g551d.extractable_allosteric_work_kb_t, uncorr.extractable_allosteric_work_kb_t
        )
        self.assertGreater(
            uncorr.extractable_allosteric_work_kb_t, w1282x.extractable_allosteric_work_kb_t
        )

        self.assertGreater(eti.mutual_information_nats / uncorr.mutual_information_nats, 4.0)
        self.assertGreater(eti.extractable_allosteric_work_kcal_mol, 1.40)

    def test_07_three_poisoned_negative_controls_fail_strictly(self) -> None:
        p_super = SagawaUedaAllostericProfile(
            name="bad_super_demon",
            description="Super-unitary Maxwell demon violation",
            i_mut_nats=2.50,
            eta_su=1.45,
            icl4_coupling_contact_fraction=0.80,
            poison_mode="super_unitary_demon_violation",
        )
        res_super = simulate_sagawa_ueda_feedback(p_super)
        self.assertFalse(res_super.is_generalized_second_law_satisfied)
        self.assertFalse(res_super.is_efficiency_within_unitary_bound)
        self.assertLess(res_super.generalized_second_law_slack_kb_t, 0.0)

        p_neg_mi = SagawaUedaAllostericProfile(
            name="bad_neg_mi",
            description="Negative mutual information violation",
            i_mut_nats=-1.60,
            eta_su=0.50,
            icl4_coupling_contact_fraction=0.50,
            poison_mode="negative_mutual_information",
        )
        res_neg_mi = simulate_sagawa_ueda_feedback(p_neg_mi)
        self.assertFalse(res_neg_mi.is_mutual_information_non_negative)
        self.assertFalse(res_neg_mi.is_generalized_second_law_satisfied)
        self.assertLess(res_neg_mi.mutual_information_nats, 0.0)

        p_broken_gj = SagawaUedaAllostericProfile(
            name="bad_gj",
            description="Broken Generalized Jarzynski equality",
            i_mut_nats=2.50,
            eta_su=0.60,
            icl4_coupling_contact_fraction=0.80,
            poison_mode="broken_generalized_jarzynski",
        )
        res_gj = simulate_sagawa_ueda_feedback(p_broken_gj)
        self.assertFalse(res_gj.is_generalized_jarzynski_valid)
        self.assertGreater(res_gj.generalized_jarzynski_error, 0.5)
        self.assertGreater(integrate_generalized_jarzynski_average(-1.0, 1.0), 2.0)


if __name__ == "__main__":
    unittest.main()
