"""Unit tests for Round 81 (T-0673): Gingrich-Horowitz First Arrival Time Bounds & Dwell-Time TUR in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.119.170601 (Gingrich TR, Horowitz JM: Fundamental Bounds on First Arrival Times for Nonequilibrium Processes, Phys. Rev. Lett. 119, 170601, 2017)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis, J. Gen. Physiol. 2006)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM, 5UAK, 8EIQ
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_gingrich_horowitz_first_passage_tur,
)


class GingrichHorowitzFirstPassageTurTests(unittest.TestCase):
    """Test suite verifying Round 81 Gingrich-Horowitz first arrival time thermodynamic bounds."""

    def setUp(self) -> None:
        self.report = evaluate_gingrich_horowitz_first_passage_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0673")
        self.assertTrue(self.report["gingrich_horowitz_first_passage_tur_bridge_verified"])
        self.assertTrue(self.report["first_arrival_bound_verified"])
        self.assertTrue(self.report["optimal_burst_timing_verified"])
        self.assertTrue(self.report["trikafta_timing_restoration_verified"])
        self.assertTrue(self.report["proof_dag_arrival_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_first_arrival_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "FIRST_ARRIVAL_TIME_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["tur_product_q_fpt"], 1.0)
        self.assertGreaterEqual(wt["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        self.assertEqual(rescued["gating_regime"], "FIRST_ARRIVAL_TIME_OPTIMAL")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["tur_product_q_fpt"], 1.0)
        self.assertGreaterEqual(rescued["uncertainty_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.85)

        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_TIMING_COLLAPSE")
        self.assertFalse(f508del["bound_satisfied"])
        self.assertLess(f508del["tur_product_q_fpt"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.20)

        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_TIMING_COLLAPSE")
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["tur_product_q_fpt"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)

        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertFalse(w1282x["bound_satisfied"])
        self.assertEqual(w1282x["tur_product_q_fpt"], 0.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["first_arrival_bound_satisfied"])
            self.assertGreaterEqual(c["tur_product_q_fpt"], 1.0)

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_gingrich_horowitz_first_passage_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
