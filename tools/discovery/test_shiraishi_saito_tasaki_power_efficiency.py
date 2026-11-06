"""Unit tests for Round 87 (T-0686): Shiraishi-Saito-Tasaki Power-Efficiency Trade-Off & CFTR NBD1-NBD2 Engine Bounds in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.117.190601 (Shiraishi N, Saito K, Tasaki H: Universal Trade-Off Relation between Power and Efficiency for Heat Engines, Phys. Rev. Lett. 117, 190601, 2016)
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
    evaluate_shiraishi_saito_tasaki_power_efficiency,
)


class ShiraishiSaitoTasakiPowerEfficiencyTests(unittest.TestCase):
    """Test suite verifying Round 87 Shiraishi-Saito-Tasaki power-efficiency trade-off bounds."""

    def setUp(self) -> None:
        self.report = evaluate_shiraishi_saito_tasaki_power_efficiency()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0686")
        self.assertTrue(self.report["shiraishi_saito_tasaki_power_efficiency_bridge_verified"])
        self.assertTrue(self.report["power_efficiency_tradeoff_bound_verified"])
        self.assertTrue(self.report["finite_time_nbd_dimerization_power_verified"])
        self.assertTrue(self.report["trikafta_power_efficiency_restoration_verified"])
        self.assertTrue(self.report["proof_dag_power_efficiency_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_power_efficiency_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "FINITE_TIME_POWER_EFFICIENCY_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertLessEqual(wt["observed_mechanical_power_kbt_s"], wt["shiraishi_power_upper_bound_kbt_s"])
        self.assertLessEqual(wt["tradeoff_saturation_ratio_r_sst"], 1.0)

        self.assertEqual(rescued["gating_regime"], "RESCUED_POWER_EFFICIENCY_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertLessEqual(rescued["observed_mechanical_power_kbt_s"], rescued["shiraishi_power_upper_bound_kbt_s"])
        self.assertLessEqual(rescued["tradeoff_saturation_ratio_r_sst"], 1.0)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["power_efficiency_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_shiraishi_saito_tasaki_power_efficiency(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
