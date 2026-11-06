"""Unit tests for Round 83 (T-0677): Garrahan Dynamical Activity First-Passage Bounds & Frenetic Gating Regularity in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevE.95.032134 (Garrahan JP: Simple bounds on fluctuations and first-passage times for counting observables, Phys. Rev. E 95, 032134, 2017)
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
    evaluate_garrahan_frenetic_first_passage_bound,
)


class GarrahanFreneticFirstPassageBoundTests(unittest.TestCase):
    """Test suite verifying Round 83 Garrahan dynamical activity first-passage bounds."""

    def setUp(self) -> None:
        self.report = evaluate_garrahan_frenetic_first_passage_bound()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0677")
        self.assertTrue(self.report["garrahan_frenetic_first_passage_bound_bridge_verified"])
        self.assertTrue(self.report["frenetic_first_passage_bound_verified"])
        self.assertTrue(self.report["optimal_frenetic_activity_verified"])
        self.assertTrue(self.report["trikafta_frenetic_restoration_verified"])
        self.assertTrue(self.report["proof_dag_frenetic_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_frenetic_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "FRENETIC_ACTIVITY_FIRST_PASSAGE_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["frenetic_product_q_garrahan"], 1.0)
        self.assertGreaterEqual(wt["relative_timing_variance_epsilon_sq"], wt["garrahan_activity_lower_bound"])

        self.assertEqual(rescued["gating_regime"], "RESCUED_FRENETIC_FIRST_PASSAGE_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["frenetic_product_q_garrahan"], 1.0)
        self.assertGreaterEqual(rescued["relative_timing_variance_epsilon_sq"], rescued["garrahan_activity_lower_bound"])

        self.assertEqual(f508del["gating_regime"], "FRENETIC_FREEZE_FIRST_PASSAGE_COLLAPSE")
        self.assertFalse(f508del["bound_satisfied"])
        self.assertLess(f508del["frenetic_product_q_garrahan"], 1.0)

        self.assertEqual(g551d["gating_regime"], "FUTILE_TRAFFIC_FIRST_PASSAGE_COLLAPSE")
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["frenetic_product_q_garrahan"], 1.0)

        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertFalse(w1282x["bound_satisfied"])
        self.assertEqual(w1282x["frenetic_product_q_garrahan"], 0.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["frenetic_bound_satisfied"])
            self.assertGreaterEqual(c["proof_frenetic_product_q"], 1.0)

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_garrahan_frenetic_first_passage_bound(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
