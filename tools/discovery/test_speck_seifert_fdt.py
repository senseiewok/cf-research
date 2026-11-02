"""Unit tests for Round 42 (T-0589): Speck-Seifert Non-Equilibrium Fluctuation-Dissipation Theorem & Velocity Field Coupling for CFTR Active Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1088/0034-4885/75/12/126001 (Seifert U: Stochastic thermodynamics, fluctuation theorems and molecular machines, Rep. Prog. Phys. 75, 126001, 2012; Sec. 3.2.1 Restoring FDT in steady states)
  - DOI:10.1103/PhysRevLett.103.090601 (Prost J, Joanny JF, Parrondo JMR: Generalized Fluctuation-Dissipation Theorem for Steady States and Applications to Molecular Motors, Phys. Rev. Lett. 103, 090601, 2009)
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
    evaluate_speck_seifert_fdt,
)


class SpeckSeifertFDTTests(unittest.TestCase):
    """Test suite verifying Round 42 Speck-Seifert non-equilibrium FDT and velocity field coupling for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_speck_seifert_fdt()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0589")
        self.assertTrue(self.report["speck_seifert_bridge_verified"])
        self.assertTrue(self.report["velocity_field_fdt_stratification_verified"])
        self.assertTrue(self.report["proof_dag_velocity_fdt_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_speck_seifert_fdt_and_gating_stratification(self) -> None:
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

        # WT CFTR: Active velocity field motor (v_ss >= 1.50 s^-1, A_FDT >= 0.65, chi_eff >= 1.50, C_v < -0.50)
        self.assertTrue(wt["generalized_fdt_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_VELOCITY_FIELD_MOTOR")
        self.assertGreaterEqual(wt["steady_state_cycling_velocity_s_inv"], 1.50)
        self.assertGreaterEqual(wt["fdt_violation_asymmetry_ratio"], 0.65)
        self.assertGreaterEqual(wt["restored_generalized_response_chi_eff"], 1.50)
        self.assertLess(wt["convective_velocity_drift_term_cv"], -0.50)

        # F508del: Passive equilibrium collapse (v_ss < 0.01 s^-1, A_FDT < 0.02, chi_eff < 0.05)
        self.assertFalse(f508del["generalized_fdt_satisfied"])
        self.assertEqual(f508del["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(f508del["steady_state_cycling_velocity_s_inv"], 0.01)
        self.assertLess(f508del["fdt_violation_asymmetry_ratio"], 0.02)
        self.assertLess(f508del["restored_generalized_response_chi_eff"], 0.05)

        # G551D: Passive equilibrium collapse (v_ss < 0.01 s^-1, A_FDT < 0.02, chi_eff < 0.05)
        self.assertFalse(g551d["generalized_fdt_satisfied"])
        self.assertEqual(g551d["gating_regime"], "PASSIVE_EQUILIBRIUM_COLLAPSE")
        self.assertLess(g551d["steady_state_cycling_velocity_s_inv"], 0.01)
        self.assertLess(g551d["fdt_violation_asymmetry_ratio"], 0.02)
        self.assertLess(g551d["restored_generalized_response_chi_eff"], 0.05)

        # Rescued CFTR (Ivacaftor / ETI): Restores active motor kinetics (v_ss >= 1.50 s^-1, A_FDT >= 0.60)
        self.assertTrue(rescued["generalized_fdt_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_VELOCITY_FIELD_MOTOR")
        self.assertGreaterEqual(rescued["steady_state_cycling_velocity_s_inv"], 1.50)
        self.assertGreaterEqual(rescued["fdt_violation_asymmetry_ratio"], 0.60)
        self.assertGreaterEqual(rescued["restored_generalized_response_chi_eff"], 1.50)

        # W1282X: Passive equilibrium APO baseline
        self.assertFalse(w1282x["generalized_fdt_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["steady_state_cycling_velocity_s_inv"], 0.0)
        self.assertEqual(w1282x["fdt_violation_asymmetry_ratio"], 0.0)
        self.assertEqual(w1282x["convective_velocity_drift_term_cv"], 0.0)
        self.assertEqual(w1282x["restored_generalized_response_chi_eff"], 0.0)

    def test_formal_proof_dag_fdt_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["fdt_proof_parity_satisfied"])
            self.assertTrue(c_prof["generalized_response_bound_satisfied"])
            self.assertGreaterEqual(c_prof["proof_inductive_response_flow"], 2.0)
            self.assertGreaterEqual(c_prof["fdt_violation_asymmetry_ratio"], 0.70)
            self.assertGreaterEqual(c_prof["deductive_implication_velocity"], 2.50)
            self.assertGreaterEqual(c_prof["topological_flux_preservation_pct"], 75.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_speck_seifert_fdt(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
