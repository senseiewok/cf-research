"""Unit tests for Round 107: Hopfield-Ninio Kinetic Proofreading, GTP Dissipation Cost, and Suppressor-tRNA Selectivity."""

import unittest
from tools.discovery.round107_hopfield_ninio_kinetic_proofreading_ptc import (
    CITATIONS,
    HYPOTHESIS,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    TITLE,
    run_falsification,
    run_falsification_checks,
    simulate_proofreading_profile,
)


class TestRound107HopfieldNinioKineticProofreading(unittest.TestCase):
    def test_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round107")
        self.assertEqual(TASK_ID, "T-0736")
        self.assertIn("Hopfield-Ninio Kinetic Proofreading", TITLE)
        self.assertIn("ACE-tRNA", HYPOTHESIS)
        self.assertGreaterEqual(len(CITATIONS), 5)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], "round107")
        self.assertTrue(ah10["passed"])

    def test_falsification_checks_all_pass(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        for name, ok in res["checks"].items():
            self.assertTrue(ok, f"Failed check: {name}")

    def test_hopfield_bound_and_dissipation(self) -> None:
        for profile in PROFILES:
            sim = simulate_proofreading_profile(profile)
            self.assertTrue(sim.is_hopfield_bound_satisfied)
            self.assertTrue(sim.is_second_law_satisfied)
            self.assertGreater(sim.gtp_dissipation_per_cycle_kbt, 0.0)
            self.assertGreaterEqual(sim.gtp_dissipation_per_peptide_kbt, sim.gtp_dissipation_per_cycle_kbt)

    def test_ace_trna_vs_aminoglycoside_selectivity(self) -> None:
        by_name = {p.name: simulate_proofreading_profile(p) for p in PROFILES}
        ace = by_name["ace_trna_trp_w1282x"]
        g418 = by_name["aminoglycoside_g418_w1282x"]
        elx = by_name["synthetic_readthrough_elx02_g542x"]
        basal = by_name["basal_near_cognate_w1282x"]

        # PTC readthrough yield
        self.assertGreater(ace.ptc_readthrough_yield_pct, 25.0)
        self.assertGreater(ace.ptc_readthrough_yield_pct, g418.ptc_readthrough_yield_pct)
        self.assertGreater(g418.ptc_readthrough_yield_pct, elx.ptc_readthrough_yield_pct)
        self.assertGreater(elx.ptc_readthrough_yield_pct, basal.ptc_readthrough_yield_pct)

        # NTC collateral misincorporation (toxic off-target)
        self.assertGreater(g418.ntc_collateral_error_pct, elx.ntc_collateral_error_pct)
        self.assertGreater(elx.ntc_collateral_error_pct, ace.ntc_collateral_error_pct)

        # Therapeutic Efficacy Index (TEI)
        self.assertGreater(ace.therapeutic_efficacy_index, 1000.0)
        self.assertGreater(ace.therapeutic_efficacy_index, elx.therapeutic_efficacy_index)
        self.assertGreater(elx.therapeutic_efficacy_index, g418.therapeutic_efficacy_index)
        self.assertEqual(basal.therapeutic_efficacy_index, 0.0)


if __name__ == "__main__":
    unittest.main()
