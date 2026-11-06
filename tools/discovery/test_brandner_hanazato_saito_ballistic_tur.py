"""Unit tests for Round 85 (T-0682): Brandner-Hanazato-Saito Universal Thermodynamic Bounds on Ballistic Transport & Channel Precision in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.120.090601 (Brandner K, Hanazato T, Saito K: Thermodynamic Bounds on Precision in Ballistic Multiterminal Transport, Phys. Rev. Lett. 120, 090601, 2018)
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
    evaluate_brandner_hanazato_saito_ballistic_tur,
)


class BrandnerHanazatoSaitoBallisticTurTests(unittest.TestCase):
    """Test suite verifying Round 85 Brandner-Hanazato-Saito ballistic multi-terminal transport TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_brandner_hanazato_saito_ballistic_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0682")
        self.assertTrue(self.report["brandner_hanazato_saito_ballistic_tur_bridge_verified"])
        self.assertTrue(self.report["ballistic_multiterminal_precision_bound_verified"])
        self.assertTrue(self.report["chloride_bicarbonate_pore_selectivity_verified"])
        self.assertTrue(self.report["trikafta_ballistic_pore_restoration_verified"])
        self.assertTrue(self.report["proof_dag_ballistic_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_ballistic_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "BALLISTIC_OPEN_PORE_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_current_relative_variance"], wt["ballistic_tur_lower_bound"])
        self.assertGreaterEqual(wt["ballistic_precision_product_q_bhs"], 2.0)

        self.assertEqual(rescued["gating_regime"], "RESCUED_BALLISTIC_PORE_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["observed_current_relative_variance"], rescued["ballistic_tur_lower_bound"])
        self.assertGreaterEqual(rescued["ballistic_precision_product_q_bhs"], 2.0)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["ballistic_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_brandner_hanazato_saito_ballistic_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
