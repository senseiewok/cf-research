"""Unit tests for Round 28 (T-0564): Schnakenberg Cycle Kinetics, Non-Equilibrium Cycle Flux Splitting & Detailed Balance Breaking in CFTR Gating.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/RevModPhys.48.571 (Schnakenberg: Network theory of microscopic and macroscopic behavior of master equation systems)
  - DOI:10.1085/jgp.200609558 (Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Cryo-EM structures of rescued CFTR)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
"""

from __future__ import annotations

import unittest
from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_schnakenberg_cycle_kinetics,
)


class SchnakenbergCycleKineticsTests(unittest.TestCase):
    """Test suite verifying Round 28 Schnakenberg cycle kinetics, cycle flux splitting, and detailed balance breaking in CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_schnakenberg_cycle_kinetics()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0564")
        self.assertTrue(self.report["schnakenberg_bridge_verified"])
        self.assertTrue(self.report["cycle_flux_splitting_verified"])
        self.assertTrue(self.report["proof_cycle_acyclicity_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_schnakenberg_cycle_kinetics_and_detailed_balance_stratification(self) -> None:
        profiles = self.report["gating_profiles"]
        self.assertIn("wt_cftr", profiles)
        self.assertIn("g551d_uncoupled_cftr", profiles)
        self.assertIn("rescued_cftr_ivacaftor", profiles)
        self.assertIn("w1282x_nonsense_baseline", profiles)

        wt = profiles["wt_cftr"]
        g551d = profiles["g551d_uncoupled_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        # WT CFTR: High cycle affinity (>15 kJ/mol), strong detailed balance breaking (>10x), active open state probability (>30%)
        self.assertTrue(wt["schnakenberg_criterion_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_DIRECTED_CYCLE")
        self.assertGreaterEqual(wt["atp_hydrolysis_affinity_kj_mol"], 18.0)
        self.assertGreaterEqual(wt["net_active_cycle_flux_s_inv"], 1.50)
        self.assertGreaterEqual(wt["detailed_balance_breaking_ratio"], 15.0)
        self.assertGreaterEqual(wt["open_probability_pct"], 30.0)

        # G551D: Uncoupled ATP hydrolysis, near-equilibrium thermal flicker, cycle flux near zero, failing Schnakenberg active criterion
        self.assertFalse(g551d["schnakenberg_criterion_satisfied"])
        self.assertEqual(g551d["gating_regime"], "PASSIVE_THERMAL_FLICKER")
        self.assertLess(g551d["atp_hydrolysis_affinity_kj_mol"], 2.0)
        self.assertLess(g551d["net_active_cycle_flux_s_inv"], 0.15)
        self.assertLess(g551d["detailed_balance_breaking_ratio"], 2.0)
        self.assertLess(g551d["open_probability_pct"], 2.0)

        # Rescued Ivacaftor: Restores non-equilibrium directed cycle kinetics, matching wild-type bounds (>85% flux recovery)
        self.assertTrue(rescued["schnakenberg_criterion_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_DIRECTED_CYCLE")
        self.assertGreaterEqual(rescued["atp_hydrolysis_affinity_kj_mol"], 18.0)
        self.assertGreaterEqual(rescued["net_active_cycle_flux_s_inv"], 1.50)
        self.assertGreaterEqual(rescued["detailed_balance_breaking_ratio"], 15.0)
        self.assertGreaterEqual(rescued["open_probability_pct"], 30.0)

        # W1282X: Passive equilibrium baseline
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["atp_hydrolysis_affinity_kj_mol"], 0.0)
        self.assertEqual(w1282x["net_active_cycle_flux_s_inv"], 0.0)
        self.assertEqual(w1282x["open_probability_pct"], 0.0)
        self.assertFalse(w1282x["schnakenberg_criterion_satisfied"])

    def test_formal_proof_cycle_acyclicity_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["schnakenberg_dag_satisfied"])
            self.assertFalse(c_prof["circular_reasoning_detected"])
            self.assertEqual(c_prof["deductive_cycle_affinity"], 0.0)
            self.assertGreater(c_prof["acyclic_reduction_flux"], 5.0)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_schnakenberg_cycle_kinetics(n_gating_regimes=2)


if __name__ == "__main__":
    unittest.main()
