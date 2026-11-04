"""Unit tests for Round 56 (T-0616): Koyuk-Seifert General Bound on Speed, Precision, and Dissipation in CFTR Biochemical Gating Cycles in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.125.260604 (Koyuk T, Seifert U: Thermodynamic Bound on the Speed of Biochemical Clocks and Oscillations, Phys. Rev. Lett. 125, 260604, 2020)
  - DOI:10.1103/PhysRevLett.122.230601 (Koyuk T, Seifert U: Operationally accessible bounds on fluctuations and dissipation in non-equilibrium steady states, Phys. Rev. Lett. 122, 230601, 2019)
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
    evaluate_koyuk_seifert_cyclic_speed_dissipation,
)


class KoyukSeifertCyclicSpeedDissipationTests(unittest.TestCase):
    """Test suite verifying Round 56 Koyuk-Seifert speed, precision, and dissipation bounds for CFTR biochemical gating cycles."""

    def setUp(self) -> None:
        self.report = evaluate_koyuk_seifert_cyclic_speed_dissipation()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0616")
        self.assertTrue(self.report["koyuk_seifert_cyclic_speed_bridge_verified"])
        self.assertTrue(self.report["cyclic_speed_bounds_verified"])
        self.assertTrue(self.report["precision_dissipation_stratification_verified"])
        self.assertTrue(self.report["proof_dag_koyuk_seifert_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_cyclic_speed_bounds_and_stratification(self) -> None:
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

        # WT CFTR: Koyuk-Seifert cyclic speed bound satisfied (tau >= lower_bound, ratio >= 1.0)
        self.assertTrue(wt["koyuk_seifert_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "CYCLIC_DISSIPATION_OPTIMAL")
        self.assertGreaterEqual(wt["cyclic_period_tau_ms"], 20.0)
        self.assertGreaterEqual(
            wt["cyclic_period_tau_ms"],
            wt["koyuk_seifert_lower_bound_ms"],
        )
        self.assertGreaterEqual(wt["cycle_precision_ratio"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores cyclic machine efficiency matching WT
        self.assertTrue(rescued["koyuk_seifert_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "CYCLIC_DISSIPATION_OPTIMAL")
        self.assertGreaterEqual(rescued["cyclic_period_tau_ms"], 20.0)
        self.assertGreaterEqual(
            rescued["cyclic_period_tau_ms"],
            rescued["koyuk_seifert_lower_bound_ms"],
        )
        self.assertGreaterEqual(rescued["cycle_precision_ratio"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Cyclic coherence collapse (tau diverges, bound violated, precision ratio < 1.0)
        self.assertFalse(f508del["koyuk_seifert_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "CYCLIC_COHERENCE_COLLAPSE")
        self.assertGreater(f508del["cyclic_period_tau_ms"], 500.0)
        self.assertLess(f508del["cycle_precision_ratio"], 1.0)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Cyclic uncoupling collapse
        self.assertFalse(g551d["koyuk_seifert_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "CYCLIC_COHERENCE_COLLAPSE")
        self.assertGreater(g551d["cyclic_period_tau_ms"], 500.0)
        self.assertLess(g551d["cycle_precision_ratio"], 1.0)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["koyuk_seifert_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["cyclic_period_tau_ms"], 0.0)
        self.assertEqual(w1282x["koyuk_seifert_lower_bound_ms"], 0.0)
        self.assertEqual(w1282x["cycle_precision_ratio"], 0.0)
        self.assertEqual(w1282x["entropy_production_rate_sigma_kb_per_ms"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_koyuk_seifert_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_koyuk_seifert_parity_satisfied"])
            self.assertGreaterEqual(c_prof["proof_cycle_duration_step"], 5.0)
            self.assertGreaterEqual(c_prof["proof_cyclic_bound"], 3.5)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_koyuk_seifert_cyclic_speed_dissipation(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
