"""Unit tests for Round 55 (T-0615): Hatano-Sasa Steady-State Thermodynamics & Dissipated Work Bounds in CFTR Gating Transitions in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.86.3463 (Hatano T, Sasa S: Steady-State Thermodynamics of Langevin Systems, Phys. Rev. Lett. 86, 3463, 2001)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
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
    evaluate_hatano_sasa_steady_state_thermodynamics,
)


class HatanoSasaSteadyStateThermodynamicsTests(unittest.TestCase):
    """Test suite verifying Round 55 Hatano-Sasa steady-state thermodynamics and dissipated work bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_hatano_sasa_steady_state_thermodynamics()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0615")
        self.assertTrue(self.report["hatano_sasa_steady_state_bridge_verified"])
        self.assertTrue(self.report["hatano_sasa_equality_verified"])
        self.assertTrue(self.report["excess_dissipated_work_stratification_verified"])
        self.assertTrue(self.report["proof_dag_hatano_sasa_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_steady_state_equality_and_dissipated_work_stratification(self) -> None:
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

        # WT CFTR: Hatano-Sasa steady-state equality satisfied (<exp(-Y)> ~ 1.0, Y >= 0)
        self.assertTrue(wt["hatano_sasa_equality_satisfied"])
        self.assertTrue(wt["hatano_sasa_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "NON_EQUILIBRIUM_STEADY_STATE_OPTIMAL")
        self.assertGreaterEqual(wt["excess_dissipated_work_y_kb_t"], 0.0)
        self.assertLessEqual(wt["excess_dissipated_work_y_kb_t"], 2.0)
        self.assertGreaterEqual(wt["steady_state_quasipotential_phi_kb_t"], 2.5)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores non-equilibrium steady state matching WT
        self.assertTrue(rescued["hatano_sasa_equality_satisfied"])
        self.assertTrue(rescued["hatano_sasa_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "NON_EQUILIBRIUM_STEADY_STATE_OPTIMAL")
        self.assertGreaterEqual(rescued["excess_dissipated_work_y_kb_t"], 0.0)
        self.assertLessEqual(rescued["excess_dissipated_work_y_kb_t"], 2.0)
        self.assertGreaterEqual(rescued["steady_state_quasipotential_phi_kb_t"], 2.5)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Severe kinetic trapping and breakdown of Hatano-Sasa transition equality
        self.assertFalse(f508del["hatano_sasa_equality_satisfied"])
        self.assertFalse(f508del["hatano_sasa_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_TRAPPING_BREAKDOWN")
        self.assertGreater(f508del["excess_dissipated_work_y_kb_t"], 4.0)
        self.assertLess(f508del["steady_state_quasipotential_phi_kb_t"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling breakdown
        self.assertFalse(g551d["hatano_sasa_equality_satisfied"])
        self.assertFalse(g551d["hatano_sasa_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "KINETIC_TRAPPING_BREAKDOWN")
        self.assertGreater(g551d["excess_dissipated_work_y_kb_t"], 4.0)
        self.assertLess(g551d["steady_state_quasipotential_phi_kb_t"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["hatano_sasa_equality_satisfied"])
        self.assertFalse(w1282x["hatano_sasa_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["excess_dissipated_work_y_kb_t"], 0.0)
        self.assertEqual(w1282x["steady_state_quasipotential_phi_kb_t"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_hatano_sasa_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_hatano_sasa_parity_satisfied"])
            self.assertGreaterEqual(c_prof["proof_excess_functional"], 1.0)
            self.assertGreaterEqual(c_prof["proof_quasipotential"], 3.5)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_hatano_sasa_steady_state_thermodynamics(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
