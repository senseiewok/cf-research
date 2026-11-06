"""Unit tests for Round 89 (T-0692): Pal-Rahav-Reuveni First-Passage Thermodynamic Uncertainty Relations under Stochastic Gating Resets in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevResearch.3.013234 (Pal A, Rahav S, Reuveni S: Thermodynamic uncertainty relation for first-passage times in reset Markov processes, Phys. Rev. Res. 3, 013234, 2021)
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
    evaluate_pal_rahav_reuveni_reset_tur,
)


class PalRahavReuveniResetTurTests(unittest.TestCase):
    """Test suite verifying Round 89 Pal-Rahav-Reuveni reset first-passage TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_pal_rahav_reuveni_reset_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0692")
        self.assertTrue(self.report["pal_rahav_reuveni_reset_tur_bridge_verified"])
        self.assertTrue(self.report["stochastic_reset_tur_bound_verified"])
        self.assertTrue(self.report["atp_unbinding_reset_rate_verified"])
        self.assertTrue(self.report["trikafta_reset_stabilization_verified"])
        self.assertTrue(self.report["proof_dag_reset_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_reset_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "STOCHASTIC_RESET_STABILIZED")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_first_passage_variance"], wt["reset_tur_lower_bound"])
        self.assertLessEqual(wt["reset_rate_r_s"], 0.1)

        self.assertEqual(rescued["gating_regime"], "RESCUED_RESET_STABILIZED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["observed_first_passage_variance"], rescued["reset_tur_lower_bound"])
        self.assertLessEqual(rescued["reset_rate_r_s"], 0.2)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["reset_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_pal_rahav_reuveni_reset_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
