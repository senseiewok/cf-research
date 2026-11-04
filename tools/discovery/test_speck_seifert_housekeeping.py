"""Unit tests for Round 60 (T-0624): Speck-Seifert Housekeeping Heat Fluctuation Theorem and Non-Equilibrium Dissipation in CFTR Gating Cycles in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0305-4470/38/34/L03 (Speck T, Seifert U: Integral fluctuation theorem for housekeeping heat, J. Phys. A: Math. Gen. 38, L581-L588, 2005)
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
    evaluate_speck_seifert_housekeeping,
)


class SpeckSeifertHousekeepingTests(unittest.TestCase):
    """Test suite verifying Round 60 Speck-Seifert housekeeping heat integral fluctuation theorem for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_speck_seifert_housekeeping()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0624")
        self.assertTrue(self.report["speck_seifert_housekeeping_bridge_verified"])
        self.assertTrue(self.report["housekeeping_ift_verified"])
        self.assertTrue(self.report["housekeeping_efficiency_stratified"])
        self.assertTrue(self.report["proof_dag_speck_seifert_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_housekeeping_ift_and_stratification(self) -> None:
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

        # WT CFTR: Speck-Seifert IFT satisfied <exp(-q_hk)> = 1.0, active housekeeping
        self.assertTrue(wt["housekeeping_ift_satisfied"])
        self.assertEqual(wt["gating_regime"], "STEADY_STATE_HOUSEKEEPING_ACTIVE")
        self.assertAlmostEqual(wt["exp_neg_q_hk_expectation"], 1.0, places=2)
        self.assertGreater(wt["housekeeping_heat_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(wt["housekeeping_coupling_efficiency"], 0.80)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores housekeeping heat
        self.assertTrue(rescued["housekeeping_ift_satisfied"])
        self.assertEqual(rescued["gating_regime"], "STEADY_STATE_HOUSEKEEPING_ACTIVE")
        self.assertAlmostEqual(rescued["exp_neg_q_hk_expectation"], 1.0, places=2)
        self.assertGreater(rescued["housekeeping_heat_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(rescued["housekeeping_coupling_efficiency"], 0.80)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Housekeeping dissipation collapse
        self.assertFalse(f508del["housekeeping_ift_satisfied"])
        self.assertEqual(f508del["gating_regime"], "HOUSEKEEPING_DISSIPATION_COLLAPSE")
        self.assertLess(f508del["housekeeping_heat_rate_kcal_mol_s"], 0.5)
        self.assertLess(f508del["housekeeping_coupling_efficiency"], 0.25)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling breakdown
        self.assertFalse(g551d["housekeeping_ift_satisfied"])
        self.assertEqual(g551d["gating_regime"], "HOUSEKEEPING_DISSIPATION_COLLAPSE")
        self.assertLess(g551d["housekeeping_heat_rate_kcal_mol_s"], 0.5)
        self.assertLess(g551d["housekeeping_coupling_efficiency"], 0.25)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["housekeeping_ift_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["housekeeping_heat_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_speck_seifert_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_speck_seifert_parity_satisfied"])
            self.assertGreater(c_prof["proof_housekeeping_cost"], 0.0)
            self.assertAlmostEqual(c_prof["proof_exp_neg_cost_expectation"], 1.0, places=2)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_speck_seifert_housekeeping(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
