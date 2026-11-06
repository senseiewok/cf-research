"""Unit tests for Round 88 (T-0691): Macieszczak-Brandner-Garrahan Unified Coherence & Thermodynamic Uncertainty Relations in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.127.220601 (Macieszczak K, Brandner K, Garrahan JP: Unified Thermodynamic Uncertainty Relations in Classical and Quantum Regimes, Phys. Rev. Lett. 127, 220601, 2021)
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
    evaluate_macieszczak_brandner_garrahan_unified_tur,
)


class MacieszczakBrandnerGarrahanUnifiedTurTests(unittest.TestCase):
    """Test suite verifying Round 88 Macieszczak-Brandner-Garrahan unified TUR bounds."""

    def setUp(self) -> None:
        self.report = evaluate_macieszczak_brandner_garrahan_unified_tur()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0691")
        self.assertTrue(self.report["macieszczak_brandner_garrahan_unified_tur_bridge_verified"])
        self.assertTrue(self.report["allosteric_coherence_bound_verified"])
        self.assertTrue(self.report["nbd_dimer_coherence_verified"])
        self.assertTrue(self.report["trikafta_coherence_restoration_verified"])
        self.assertTrue(self.report["proof_dag_coherence_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_coherence_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "ALLOSTERIC_COHERENT_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertGreaterEqual(wt["observed_relative_variance"], wt["unified_tur_lower_bound"])
        self.assertGreaterEqual(wt["coherence_factor_c"], 0.8)

        self.assertEqual(rescued["gating_regime"], "RESCUED_COHERENCE_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertGreaterEqual(rescued["observed_relative_variance"], rescued["unified_tur_lower_bound"])
        self.assertGreaterEqual(rescued["coherence_factor_c"], 0.7)

        self.assertFalse(f508del["bound_satisfied"])
        self.assertFalse(g551d["bound_satisfied"])
        self.assertFalse(w1282x["bound_satisfied"])

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            self.assertTrue(challenges[c_id]["parity_verified"])
            self.assertTrue(challenges[c_id]["coherence_bound_satisfied"])

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_macieszczak_brandner_garrahan_unified_tur(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
