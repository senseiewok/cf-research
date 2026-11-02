"""Unit tests for Round 43 (T-0591): Kawai-Parrondo-Van den Broeck (KPV) Dissipation-Information Relative Entropy Equality for CFTR Trajectories in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012; Sec. 2.4 Relative entropy and dissipated work)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_kpv_relative_entropy,
)


class KPVRelativeEntropyTests(unittest.TestCase):
    """Test suite verifying Round 43 Kawai-Parrondo-Van den Broeck relative entropy equality for CFTR gating trajectories."""

    def setUp(self) -> None:
        self.report = evaluate_kpv_relative_entropy()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0591")
        self.assertTrue(self.report["kpv_bridge_verified"])
        self.assertTrue(self.report["relative_entropy_dissipation_stratification_verified"])
        self.assertTrue(self.report["proof_dag_relative_entropy_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_kpv_relative_entropy_and_gating_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("f508del_uncorrected", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: Active irreversible gating cycle (D(P || P_tilde) >= 4.5 kB, <W_diss> >= 18.0 kBT, ratio >= 50.0)
        self.assertTrue(wt["kpv_equality_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_IRREVERSIBLE_CYCLE")
        self.assertGreaterEqual(wt["relative_entropy_divergence_kb"], 4.5)
        self.assertGreaterEqual(wt["dissipated_work_kbt"], 18.0)
        self.assertGreater(wt["total_entropy_production_rate_w_k"], 1e-20)
        self.assertGreater(wt["irreversibility_ratio_forward_backward"], 50.0)

        # F508del: Passive equilibrium collapse (D(P || P_tilde) < 0.20 kB, <W_diss> < 1.0 kBT)
        self.assertFalse(f508del["kpv_equality_satisfied"])
        self.assertEqual(f508del["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(f508del["relative_entropy_divergence_kb"], 0.20)
        self.assertLess(f508del["dissipated_work_kbt"], 1.0)
        self.assertLess(f508del["total_entropy_production_rate_w_k"], 1e-21)
        self.assertLess(f508del["irreversibility_ratio_forward_backward"], 2.0)

        # G551D: Degenerate futile dissipation trap (D(P || P_tilde) < 0.20 kB, <W_diss> < 1.0 kBT)
        self.assertFalse(g551d["kpv_equality_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_FUTILE_TRAP")
        self.assertLess(g551d["relative_entropy_divergence_kb"], 0.20)
        self.assertLess(g551d["dissipated_work_kbt"], 1.0)
        self.assertLess(g551d["total_entropy_production_rate_w_k"], 1e-21)
        self.assertLess(g551d["irreversibility_ratio_forward_backward"], 2.0)

        # Rescued CFTR (Ivacaftor / ETI): Restores active irreversible gating (D(P || P_tilde) >= 4.0 kB, <W_diss> >= 16.0 kBT)
        self.assertTrue(rescued["kpv_equality_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_IRREVERSIBLE_CYCLE")
        self.assertGreaterEqual(rescued["relative_entropy_divergence_kb"], 4.0)
        self.assertGreaterEqual(rescued["dissipated_work_kbt"], 16.0)
        self.assertGreater(rescued["total_entropy_production_rate_w_k"], 1e-20)
        self.assertGreater(rescued["irreversibility_ratio_forward_backward"], 50.0)

        # W1282X: Passive equilibrium APO baseline (D(P || P_tilde) == 0, <W_diss> == 0, detailed balance preserved)
        self.assertFalse(w1282x["kpv_equality_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["relative_entropy_divergence_kb"], 0.0)
        self.assertEqual(w1282x["dissipated_work_kbt"], 0.0)
        self.assertEqual(w1282x["total_entropy_production_rate_w_k"], 0.0)
        self.assertEqual(w1282x["irreversibility_ratio_forward_backward"], 1.0)

    def test_formal_proof_dag_kpv_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["kpv_proof_parity_satisfied"])
            self.assertTrue(c_prof["inductive_information_bound_satisfied"])
            self.assertGreaterEqual(c_prof["proof_information_gain_nats"], 3.5)
            self.assertGreaterEqual(c_prof["proof_relative_entropy_divergence"], 3.5)
            self.assertGreaterEqual(c_prof["deductive_arrow_of_time_asymmetry"], 0.70)
            self.assertGreaterEqual(c_prof["topological_flux_preservation_pct"], 75.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_kpv_relative_entropy(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
