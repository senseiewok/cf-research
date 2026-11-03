"""Unit tests for Round 47 (T-0599): Hartich-Barato-Seifert Bipartite Information Thermodynamics & Subsystem Entropy Bounds in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/1742-5468/2014/10/P10016 (Hartich D, Barato AC, Seifert U: Stochastic thermodynamics of bipartite systems: transfer entropy and the second law, J. Stat. Mech. 2014, P10016)
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012)
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
    evaluate_hartich_barato_seifert_bipartite,
)


class HartichBaratoSeifertBipartiteTests(unittest.TestCase):
    """Test suite verifying Round 47 Hartich-Barato-Seifert bipartite information thermodynamics and subsystem entropy bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_hartich_barato_seifert_bipartite()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0599")
        self.assertTrue(self.report["hartich_barato_seifert_bipartite_bridge_verified"])
        self.assertTrue(self.report["subsystem_entropy_bounds_verified"])
        self.assertTrue(self.report["information_flow_stratification_verified"])
        self.assertTrue(self.report["proof_dag_bipartite_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_bipartite_entropy_bounds_and_information_flow_stratification(self) -> None:
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

        # WT CFTR: Coherent bipartite information flow (I_X->Y >= 2.5 kB/s, sigma_tot >= 4.0 kB/s)
        self.assertTrue(wt["bipartite_second_law_satisfied"])
        self.assertEqual(wt["gating_regime"], "COHERENT_BIPARTITE_TRANSDUCTION")
        self.assertGreaterEqual(wt["information_flow_x_to_y_kb_per_s"], 2.5)
        self.assertGreaterEqual(wt["total_dissipation_rate_sigma_tot_kb_per_s"], 4.0)
        self.assertGreaterEqual(wt["subsystem_x_entropy_rate_kb_per_s"], 3.0)
        self.assertGreaterEqual(wt["subsystem_y_entropy_rate_kb_per_s"], 1.0)
        self.assertGreaterEqual(wt["nbd_tmd_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores bipartite transfer entropy and allosteric information flow matching WT
        self.assertTrue(rescued["bipartite_second_law_satisfied"])
        self.assertEqual(rescued["gating_regime"], "COHERENT_BIPARTITE_TRANSDUCTION")
        self.assertGreaterEqual(rescued["information_flow_x_to_y_kb_per_s"], 2.5)
        self.assertGreaterEqual(rescued["total_dissipation_rate_sigma_tot_kb_per_s"], 4.0)
        self.assertGreaterEqual(rescued["subsystem_x_entropy_rate_kb_per_s"], 3.0)
        self.assertGreaterEqual(rescued["subsystem_y_entropy_rate_kb_per_s"], 1.0)
        self.assertGreaterEqual(rescued["nbd_tmd_coupling_fidelity"], 0.90)

        # F508del: Broken bipartite transmission (I_X->Y -> 0, futile uncoupled dissipation without pore opening)
        self.assertFalse(f508del["bipartite_second_law_satisfied"])
        self.assertEqual(f508del["gating_regime"], "UNCOUPLED_FUTILE_DISSIPATION")
        self.assertLess(f508del["information_flow_x_to_y_kb_per_s"], 0.10)
        self.assertLess(f508del["total_dissipation_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(f508del["nbd_tmd_coupling_fidelity"], 0.10)

        # G551D: Broken bipartite transmission (I_X->Y -> 0, uncoupled ATP turnover in NBD subsystem without gating)
        self.assertFalse(g551d["bipartite_second_law_satisfied"])
        self.assertEqual(g551d["gating_regime"], "UNCOUPLED_FUTILE_DISSIPATION")
        self.assertLess(g551d["information_flow_x_to_y_kb_per_s"], 0.10)
        self.assertLess(g551d["total_dissipation_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(g551d["nbd_tmd_coupling_fidelity"], 0.10)

        # W1282X: Passive equilibrium baseline (zero information flow, zero subsystem dissipation)
        self.assertFalse(w1282x["bipartite_second_law_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["information_flow_x_to_y_kb_per_s"], 0.0)
        self.assertEqual(w1282x["information_flow_y_to_x_kb_per_s"], 0.0)
        self.assertEqual(w1282x["subsystem_x_entropy_rate_kb_per_s"], 0.0)
        self.assertEqual(w1282x["subsystem_y_entropy_rate_kb_per_s"], 0.0)
        self.assertEqual(w1282x["total_dissipation_rate_sigma_tot_kb_per_s"], 0.0)
        self.assertEqual(w1282x["nbd_tmd_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_bipartite_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_bipartite_decomposition_satisfied"])
            self.assertTrue(c_prof["hartich_barato_seifert_parity_satisfied"])
            self.assertGreaterEqual(c_prof["mutual_information_flow_bits_per_step"], 4.0)
            self.assertGreaterEqual(c_prof["subsystem_x_modularity_lifting_rate"], 4.0)
            self.assertGreaterEqual(c_prof["subsystem_y_level_lowering_rate"], 3.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_hartich_barato_seifert_bipartite(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
