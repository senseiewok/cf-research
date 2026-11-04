"""Unit tests for Round 59 (T-0623): Lebowitz-Spohn Stochastic Action Functional and Gallavotti-Cohen Symmetry in CFTR Gating Transitions in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1023/A:1004589714161 (Lebowitz JL, Spohn H: A Gallavotti-Cohen-type symmetry in the large deviation functional for stochastic dynamics, J. Stat. Phys. 95, 333-365, 1999)
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
    evaluate_lebowitz_spohn_action,
)


class LebowitzSpohnActionTests(unittest.TestCase):
    """Test suite verifying Round 59 Lebowitz-Spohn action functional and Gallavotti-Cohen symmetry for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_lebowitz_spohn_action()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0623")
        self.assertTrue(self.report["lebowitz_spohn_action_bridge_verified"])
        self.assertTrue(self.report["gallavotti_cohen_symmetry_verified"])
        self.assertTrue(self.report["action_functional_stratification_verified"])
        self.assertTrue(self.report["proof_dag_lebowitz_spohn_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gallavotti_cohen_symmetry_and_stratification(self) -> None:
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

        # WT CFTR: Gallavotti-Cohen symmetry slope = 1.0, active cycle
        self.assertTrue(wt["gallavotti_cohen_symmetry_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_GALLAVOTTI_COHEN_CYCLE")
        self.assertAlmostEqual(wt["gallavotti_cohen_slope"], 1.0, places=2)
        self.assertGreater(wt["trajectory_action_rate_s_inv"], 4.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores Gallavotti-Cohen symmetry
        self.assertTrue(rescued["gallavotti_cohen_symmetry_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_GALLAVOTTI_COHEN_CYCLE")
        self.assertAlmostEqual(rescued["gallavotti_cohen_slope"], 1.0, places=2)
        self.assertGreater(rescued["trajectory_action_rate_s_inv"], 4.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Severe equilibrium flicker collapse
        self.assertFalse(f508del["gallavotti_cohen_symmetry_satisfied"])
        self.assertEqual(f508del["gating_regime"], "EQUILIBRIUM_FLICKER_COLLAPSE")
        self.assertLess(f508del["trajectory_action_rate_s_inv"], 0.5)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling breakdown
        self.assertFalse(g551d["gallavotti_cohen_symmetry_satisfied"])
        self.assertEqual(g551d["gating_regime"], "EQUILIBRIUM_FLICKER_COLLAPSE")
        self.assertLess(g551d["trajectory_action_rate_s_inv"], 0.5)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["gallavotti_cohen_symmetry_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["trajectory_action_rate_s_inv"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_lebowitz_spohn_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_lebowitz_spohn_parity_satisfied"])
            self.assertGreater(c_prof["proof_directional_action"], 0.0)
            self.assertAlmostEqual(c_prof["proof_symmetry_slope"], 1.0, places=2)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_lebowitz_spohn_action(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
