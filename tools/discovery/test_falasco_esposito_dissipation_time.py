"""Unit tests for Round 57 (T-0619): Falasco-Esposito Dissipation-Time Uncertainty Relation in CFTR Non-Equilibrium Transitions in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.125.120604 (Falasco G, Esposito M: Dissipation-Time Uncertainty Relation for Driven Systems, Phys. Rev. Lett. 125, 120604, 2020)
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
    evaluate_falasco_esposito_dissipation_time,
)


class FalascoEspositoDissipationTimeTests(unittest.TestCase):
    """Test suite verifying Round 57 Falasco-Esposito dissipation-time uncertainty relation for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_falasco_esposito_dissipation_time()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0619")
        self.assertTrue(self.report["falasco_esposito_dissipation_time_bridge_verified"])
        self.assertTrue(self.report["dissipation_time_bounds_verified"])
        self.assertTrue(self.report["stratification_verified"])
        self.assertTrue(self.report["proof_dag_falasco_esposito_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_dissipation_time_bounds_and_stratification(self) -> None:
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

        # WT CFTR: Dissipation-time uncertainty relation satisfied (sigma * tau >= lower_bound, ratio > 1.0)
        self.assertTrue(wt["falasco_esposito_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "DISSIPATION_TIME_OPTIMAL")
        self.assertGreater(wt["dissipation_time_product"], wt["falasco_esposito_lower_bound"])
        self.assertGreaterEqual(wt["bound_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores dissipation-time optimality
        self.assertTrue(rescued["falasco_esposito_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "DISSIPATION_TIME_OPTIMAL")
        self.assertGreater(rescued["dissipation_time_product"], rescued["falasco_esposito_lower_bound"])
        self.assertGreaterEqual(rescued["bound_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Severe kinetic trapping violates trade-off bound (ratio < 1.0)
        self.assertFalse(f508del["falasco_esposito_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_TRAPPING_BREAKDOWN")
        self.assertLess(f508del["bound_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling breakdown
        self.assertFalse(g551d["falasco_esposito_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "KINETIC_TRAPPING_BREAKDOWN")
        self.assertLess(g551d["bound_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["falasco_esposito_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["dissipation_time_product"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_falasco_esposito_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_falasco_esposito_parity_satisfied"])
            self.assertGreater(c_prof["proof_transition_steps"], 0.0)
            self.assertGreater(c_prof["proof_dissipation_bound"], 0.0)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_falasco_esposito_dissipation_time(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
