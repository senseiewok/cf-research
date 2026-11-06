"""Unit tests for Round 86 (T-0685): Proesmans-Van den Broeck Discrete-Time TUR & Stroboscopic CFTR Gating Bounds in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1209/0295-5075/119/20001 (Proesmans K, Van den Broeck C: Discrete-time thermodynamic uncertainty relation, EPL 119, 20001, 2017)
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
    evaluate_proesmans_discrete_time_tur,
)


class ProesmansDiscreteTimeTurTests(unittest.TestCase):
    """Test suite verifying Round 86 Proesmans-Van den Broeck discrete-time TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_proesmans_discrete_time_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0685")
        self.assertTrue(self.report["proesmans_discrete_time_tur_bridge_verified"])
        self.assertTrue(self.report["discrete_exponential_fluctuation_bound_verified"])
        self.assertTrue(self.report["stroboscopic_gating_constancy_verified"])
        self.assertTrue(self.report["trikafta_discrete_bin_restoration_verified"])
        self.assertTrue(self.report["proof_dag_discrete_step_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_discrete_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "STROBOSCOPIC_DISCRETE_STEP_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_step_relative_variance"], wt["discrete_tur_lower_bound"])
        self.assertGreaterEqual(wt["discrete_constancy_product_q_disc"], 2.0)

        self.assertEqual(rescued["gating_regime"], "RESCUED_DISCRETE_GATING_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["observed_step_relative_variance"], rescued["discrete_tur_lower_bound"])
        self.assertGreaterEqual(rescued["discrete_constancy_product_q_disc"], 2.0)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["discrete_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_proesmans_discrete_time_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
