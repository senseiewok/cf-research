"""Unit tests for Round 84 (T-0681): Koyuk-Seifert Periodically Driven Thermodynamic Uncertainty Relations in Cyclic CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.122.230601 (Koyuk T, Seifert U: Operationally Accessible Bound on Fluctuations and Entropy Production in Periodically Driven Systems, Phys. Rev. Lett. 122, 230601, 2019)
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
    evaluate_koyuk_seifert_periodic_tur,
)


class KoyukSeifertPeriodicTurTests(unittest.TestCase):
    """Test suite verifying Round 84 Koyuk-Seifert periodically driven TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_koyuk_seifert_periodic_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0681")
        self.assertTrue(self.report["koyuk_seifert_periodic_tur_bridge_verified"])
        self.assertTrue(self.report["periodic_current_fluctuation_bound_verified"])
        self.assertTrue(self.report["cyclic_atp_driving_constancy_verified"])
        self.assertTrue(self.report["trikafta_periodic_restoration_verified"])
        self.assertTrue(self.report["proof_dag_periodic_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_periodic_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "PERIODIC_ATP_CYCLE_COHERENT")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_relative_variance"], wt["periodic_tur_lower_bound"])
        self.assertGreaterEqual(wt["periodic_constancy_product_q_tau"], 2.0)

        self.assertEqual(rescued["gating_regime"], "RESCUED_PERIODIC_CYCLE_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["observed_relative_variance"], rescued["periodic_tur_lower_bound"])
        self.assertGreaterEqual(rescued["periodic_constancy_product_q_tau"], 2.0)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["periodic_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_koyuk_seifert_periodic_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
