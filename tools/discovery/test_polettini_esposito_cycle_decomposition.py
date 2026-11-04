"""Unit tests for Round 61 (T-0627): Polettini-Esposito Non-Convex Relative Entropy and Cycle Dissipation Decomposition in CFTR Gating Cycles in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevE.88.012112 (Polettini M, Esposito M: Nonconvexity of the relative entropy and geometric decomposition of entropy production, Phys. Rev. E 88, 012112, 2013)
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
    evaluate_polettini_esposito_cycle_decomposition,
)


class PolettiniEspositoCycleDecompositionTests(unittest.TestCase):
    """Test suite verifying Round 61 Polettini-Esposito cycle dissipation decomposition for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_polettini_esposito_cycle_decomposition()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0627")
        self.assertTrue(self.report["polettini_esposito_cycle_bridge_verified"])
        self.assertTrue(self.report["cycle_decomposition_verified"])
        self.assertTrue(self.report["cycle_affinity_stratified"])
        self.assertTrue(self.report["proof_dag_cycle_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_cycle_decomposition_and_stratification(self) -> None:
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

        # WT CFTR: Concentrated cycle dissipation, high planar affinity >= 18 kJ/mol
        self.assertTrue(wt["cycle_decomposition_satisfied"])
        self.assertEqual(wt["gating_regime"], "PLANAR_CYCLE_DISSIPATION_CONCENTRATED")
        self.assertGreaterEqual(wt["cycle_affinity_kj_mol"], 18.0)
        self.assertGreater(wt["cycle_dissipation_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(wt["cycle_dissipation_ratio"], 0.85)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores cycle affinity and concentrated dissipation
        self.assertTrue(rescued["cycle_decomposition_satisfied"])
        self.assertEqual(rescued["gating_regime"], "PLANAR_CYCLE_DISSIPATION_CONCENTRATED")
        self.assertGreaterEqual(rescued["cycle_affinity_kj_mol"], 18.0)
        self.assertGreater(rescued["cycle_dissipation_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(rescued["cycle_dissipation_ratio"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Cycle affinity collapse
        self.assertFalse(f508del["cycle_decomposition_satisfied"])
        self.assertEqual(f508del["gating_regime"], "CYCLE_AFFINITY_COLLAPSE")
        self.assertLess(f508del["cycle_affinity_kj_mol"], 5.0)
        self.assertLess(f508del["cycle_dissipation_rate_kcal_mol_s"], 0.5)
        self.assertLess(f508del["cycle_dissipation_ratio"], 0.25)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic cycle breakdown
        self.assertFalse(g551d["cycle_decomposition_satisfied"])
        self.assertEqual(g551d["gating_regime"], "CYCLE_AFFINITY_COLLAPSE")
        self.assertLess(g551d["cycle_affinity_kj_mol"], 5.0)
        self.assertLess(g551d["cycle_dissipation_rate_kcal_mol_s"], 0.5)
        self.assertLess(g551d["cycle_dissipation_ratio"], 0.25)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["cycle_decomposition_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["cycle_affinity_kj_mol"], 0.0)
        self.assertEqual(w1282x["cycle_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_cycle_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_cycle_parity_satisfied"])
            self.assertGreater(c_prof["proof_fundamental_cycle_affinity"], 10.0)
            self.assertGreaterEqual(c_prof["proof_cycle_dissipation_ratio"], 0.80)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_polettini_esposito_cycle_decomposition(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
