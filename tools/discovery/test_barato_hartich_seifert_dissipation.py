"""Unit tests for Round 50 (T-0605): Barato-Hartich-Seifert Information Thermodynamics & Dissipation Bounds in CFTR Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1742-5468/2014/10/P10016 (Hartich D, Barato AC, Seifert U: Stochastic thermodynamics of bipartite systems: transfer entropy and the second law, J. Stat. Mech. P10016, 2014)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of Δ508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_barato_hartich_seifert_dissipation,
)


class BaratoHartichSeifertDissipationTests(unittest.TestCase):
    """Test suite verifying Round 50 Barato-Hartich-Seifert information thermodynamics and dissipation bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_barato_hartich_seifert_dissipation()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0605")
        self.assertTrue(self.report["barato_hartich_seifert_dissipation_bridge_verified"])
        self.assertTrue(self.report["bipartite_dissipation_bounds_verified"])
        self.assertTrue(self.report["information_flow_stratification_verified"])
        self.assertTrue(self.report["proof_dag_dissipation_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_bipartite_dissipation_bounds_and_information_flow_stratification(self) -> None:
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

        # WT CFTR: Barato-Hartich-Seifert dissipation bounds satisfied (sigma_X >= I_X->Y, sigma_Y >= -I_X->Y)
        # Information flow I_X->Y >= 1.5 kB/s, coherent bipartite transduction
        self.assertTrue(wt["bhs_dissipation_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "COHERENT_BIPARTITE_DISSIPATION")
        self.assertGreaterEqual(wt["information_flow_rate_i_x_to_y_kb_per_s"], 1.5)
        self.assertGreaterEqual(
            wt["subsystem_x_dissipation_rate_sigma_x_kb_per_s"],
            wt["information_flow_rate_i_x_to_y_kb_per_s"],
        )
        self.assertGreaterEqual(
            wt["subsystem_y_dissipation_rate_sigma_y_kb_per_s"],
            -wt["information_flow_rate_i_x_to_y_kb_per_s"],
        )
        self.assertGreaterEqual(wt["total_dissipation_rate_sigma_tot_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores positive information flow and dissipation bounds matching WT
        self.assertTrue(rescued["bhs_dissipation_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "COHERENT_BIPARTITE_DISSIPATION")
        self.assertGreaterEqual(rescued["information_flow_rate_i_x_to_y_kb_per_s"], 1.5)
        self.assertGreaterEqual(
            rescued["subsystem_x_dissipation_rate_sigma_x_kb_per_s"],
            rescued["information_flow_rate_i_x_to_y_kb_per_s"],
        )
        self.assertGreaterEqual(
            rescued["subsystem_y_dissipation_rate_sigma_y_kb_per_s"],
            -rescued["information_flow_rate_i_x_to_y_kb_per_s"],
        )
        self.assertGreaterEqual(rescued["total_dissipation_rate_sigma_tot_kb_per_s"], 4.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Degenerate information decoupling (I_X->Y -> 0, uncoupled dissipation failure)
        self.assertFalse(f508del["bhs_dissipation_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "DEGENERATE_INFORMATION_DECOUPLING")
        self.assertLess(f508del["information_flow_rate_i_x_to_y_kb_per_s"], 0.10)
        self.assertLess(f508del["total_dissipation_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.10)

        # G551D: Degenerate information decoupling (I_X->Y -> 0, futile uncoupled catalytic cycle)
        self.assertFalse(g551d["bhs_dissipation_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "DEGENERATE_INFORMATION_DECOUPLING")
        self.assertLess(g551d["information_flow_rate_i_x_to_y_kb_per_s"], 0.10)
        self.assertLess(g551d["total_dissipation_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)

        # W1282X: Passive equilibrium baseline (sigma_X = sigma_Y = I = 0)
        self.assertFalse(w1282x["bhs_dissipation_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["subsystem_x_dissipation_rate_sigma_x_kb_per_s"], 0.0)
        self.assertEqual(w1282x["subsystem_y_dissipation_rate_sigma_y_kb_per_s"], 0.0)
        self.assertEqual(w1282x["total_dissipation_rate_sigma_tot_kb_per_s"], 0.0)
        self.assertEqual(w1282x["information_flow_rate_i_x_to_y_kb_per_s"], 0.0)
        self.assertEqual(w1282x["information_flow_rate_i_y_to_x_kb_per_s"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_dissipation_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_bipartite_dissipation_satisfied"])
            self.assertTrue(c_prof["bhs_parity_satisfied"])
            self.assertGreaterEqual(c_prof["subsystem_x_modularity_dissipation_rate"], 4.0)
            self.assertGreaterEqual(c_prof["subsystem_y_ribet_lowering_rate"], 3.0)
            self.assertGreaterEqual(c_prof["mutual_information_flow_bits_per_step"], 4.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_barato_hartich_seifert_dissipation(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
