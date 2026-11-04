"""Unit tests for Round 58 (T-0620): Dechant-Sasa Fluctuation-Response Inequality and Out-of-Equilibrium Susceptibility in CFTR Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1073/pnas.1918386117 (Dechant A, Sasa S: Fluctuation-response inequality out of equilibrium, Proc. Natl. Acad. Sci. U.S.A. 117, 6430-6436, 2020)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of Δ508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_dechant_sasa_fluctuation_response,
)


class DechantSasaFluctuationResponseTests(unittest.TestCase):
    """Test suite verifying Round 58 Dechant-Sasa fluctuation-response inequality for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_dechant_sasa_fluctuation_response()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0620")
        self.assertTrue(self.report["dechant_sasa_fluctuation_response_bridge_verified"])
        self.assertTrue(self.report["fluctuation_response_inequality_verified"])
        self.assertTrue(self.report["susceptibility_stratification_verified"])
        self.assertTrue(self.report["proof_dag_dechant_sasa_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_fluctuation_response_bounds_and_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("f508del_uncorrected", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Fluctuation-response inequality satisfied (R^2 <= 0.5 * Var * sigma, ratio <= 1.0)
        self.assertTrue(wt["dechant_sasa_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "RESPONSE_DISSIPATION_OPTIMAL")
        self.assertLessEqual(wt["squared_response"], wt["response_upper_bound"])
        self.assertLessEqual(wt["fluctuation_response_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores response-dissipation optimality
        self.assertTrue(rescued["dechant_sasa_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "RESPONSE_DISSIPATION_OPTIMAL")
        self.assertLessEqual(rescued["squared_response"], rescued["response_upper_bound"])
        self.assertLessEqual(rescued["fluctuation_response_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Severe kinetic trapping leads to susceptibility collapse (ratio > 1.0)
        self.assertFalse(f508del["dechant_sasa_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "SUSCEPTIBILITY_COLLAPSE")
        self.assertGreater(f508del["fluctuation_response_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling susceptibility collapse
        self.assertFalse(g551d["dechant_sasa_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "SUSCEPTIBILITY_COLLAPSE")
        self.assertGreater(g551d["fluctuation_response_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["dechant_sasa_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["allosteric_susceptibility_r"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_dechant_sasa_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_dechant_sasa_parity_satisfied"])
            self.assertGreater(c_prof["proof_susceptibility"], 0.0)
            self.assertGreater(c_prof["proof_response_bound"], 0.0)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_dechant_sasa_fluctuation_response(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
