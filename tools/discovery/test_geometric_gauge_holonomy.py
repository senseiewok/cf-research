"""Unit tests for Round 37 (T-0581): Geometrical Thermodynamic Gauge Theory & Berry Holonomy for CFTR Catalytic Cycles in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1063/1.2738753 (Maes & Netočný: Minimum entropy production principle from a dynamical fluctuation law, J. Math. Phys. 48, 053306, 2007)
  - DOI:10.1085/jgp.200609558 (Csanády et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Cryo-EM structures of rescued human CFTR, PDB 8EIQ)
  - DOI:10.1103/PhysRevLett.86.3463 (Hatano & Sasa: Steady-State Thermodynamics of Transitions between Non-Equilibrium Steady States)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import math
import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_geometric_gauge_holonomy,
)


class GeometricGaugeHolonomyTests(unittest.TestCase):
    """Test suite verifying Round 37 Geometrical Thermodynamic Gauge Theory & Berry Holonomy in CFTR catalytic cycles."""

    def setUp(self) -> None:
        self.report = evaluate_geometric_gauge_holonomy()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0581")
        self.assertTrue(self.report["geometric_gauge_bridge_verified"])
        self.assertTrue(self.report["gauge_curvature_stratification_verified"])
        self.assertTrue(self.report["proof_dag_holonomy_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gauge_curvature_and_holonomy_stratification(self) -> None:
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

        # WT CFTR: Non-zero quantized geometric phase (gamma >= pi/2), directional anion pumping
        self.assertTrue(wt["geometric_phase_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "QUANTIZED_BERRY_HOLONOMY")
        self.assertGreaterEqual(wt["berry_phase_gamma_rad"], math.pi / 2.0)
        self.assertGreaterEqual(wt["gauge_curvature_f_norm"], 2.0)
        self.assertGreaterEqual(wt["holonomic_pumping_efficiency_pct"], 70.0)
        self.assertGreaterEqual(wt["anion_flux_directionality_ratio"], 0.80)
        self.assertGreaterEqual(wt["opening_transition_rate_k_s_inv"], 1.50)

        # F508del: Trivial gauge curvature (gamma ~ 0), broken cyclic transport
        self.assertFalse(f508del["geometric_phase_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "TRIVIAL_GAUGE_COLLAPSE")
        self.assertLess(f508del["berry_phase_gamma_rad"], 0.10)
        self.assertLess(f508del["gauge_curvature_f_norm"], 0.10)
        self.assertLess(f508del["holonomic_pumping_efficiency_pct"], 10.0)
        self.assertLess(f508del["anion_flux_directionality_ratio"], 0.15)
        self.assertLess(f508del["opening_transition_rate_k_s_inv"], 0.01)

        # G551D: Trivial gauge curvature (gamma ~ 0), broken cyclic transport
        self.assertFalse(g551d["geometric_phase_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "TRIVIAL_GAUGE_COLLAPSE")
        self.assertLess(g551d["berry_phase_gamma_rad"], 0.10)
        self.assertLess(g551d["gauge_curvature_f_norm"], 0.10)
        self.assertLess(g551d["holonomic_pumping_efficiency_pct"], 10.0)
        self.assertLess(g551d["anion_flux_directionality_ratio"], 0.15)
        self.assertLess(g551d["opening_transition_rate_k_s_inv"], 0.01)

        # Rescued CFTR (Ivacaftor / ETI): Restores non-trivial Berry holonomy matching WT gating cycles
        self.assertTrue(rescued["geometric_phase_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "QUANTIZED_BERRY_HOLONOMY")
        self.assertGreaterEqual(rescued["berry_phase_gamma_rad"], math.pi / 2.0)
        self.assertGreaterEqual(rescued["gauge_curvature_f_norm"], 2.0)
        self.assertGreaterEqual(rescued["holonomic_pumping_efficiency_pct"], 60.0)
        self.assertGreaterEqual(rescued["anion_flux_directionality_ratio"], 0.75)
        self.assertGreaterEqual(rescued["opening_transition_rate_k_s_inv"], 1.50)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["geometric_phase_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["berry_phase_gamma_rad"], 0.0)
        self.assertEqual(w1282x["gauge_curvature_f_norm"], 0.0)
        self.assertEqual(w1282x["opening_transition_rate_k_s_inv"], 0.0)

    def test_formal_proof_dag_holonomy_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["geometric_gauge_proof_parity_satisfied"])
            self.assertTrue(c_prof["holonomy_bound_satisfied"])
            self.assertGreaterEqual(c_prof["proof_dag_holonomy_gamma"], math.pi / 2.0)
            self.assertEqual(c_prof["topological_invariant_winding"], 1)
            self.assertGreaterEqual(c_prof["derivation_barrier_reduction_pct"], 50.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_geometric_gauge_holonomy(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
