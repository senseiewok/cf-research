"""Unit tests for Round 54 (T-0612): Esposito-Van den Broeck Three-Way Entropy Production Decomposition & Housekeeping Dissipation in CFTR Gating in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.104.090601 (Esposito M, Van den Broeck C: Three faces of the second law. I. Master equations, Phys. Rev. Lett. 104, 090601, 2010)
  - DOI:10.1103/PhysRevE.82.011143 (Esposito M, Van den Broeck C: Three faces of the second law. II. Fokker-Planck formulation, Phys. Rev. E 82, 011143, 2010)
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
    evaluate_esposito_van_den_broeck_entropy,
)


class EspositoVanDenBroeckEntropyTests(unittest.TestCase):
    """Test suite verifying Round 54 Esposito-Van den Broeck three-way entropy production decomposition and housekeeping dissipation for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_esposito_van_den_broeck_entropy()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0612")
        self.assertTrue(self.report["esposito_van_den_broeck_entropy_bridge_verified"])
        self.assertTrue(self.report["three_way_decomposition_verified"])
        self.assertTrue(self.report["housekeeping_dissipation_stratification_verified"])
        self.assertTrue(self.report["proof_dag_entropy_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_three_way_decomposition_and_housekeeping_stratification(self) -> None:
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

        # WT CFTR: Three-way decomposition satisfied (sigma_tot = sigma_hk + sigma_ex, sigma_hk >= 1.5, housekeeping fraction >= 0.50)
        self.assertTrue(wt["decomposition_equality_satisfied"])
        self.assertEqual(wt["gating_regime"], "ACTIVE_STEADY_STATE_DISSIPATION")
        self.assertGreaterEqual(wt["total_entropy_production_rate_sigma_tot_kb_per_s"], 2.5)
        self.assertGreaterEqual(wt["housekeeping_entropy_rate_sigma_hk_kb_per_s"], 1.5)
        self.assertGreaterEqual(wt["excess_entropy_rate_sigma_ex_kb_per_s"], 0.8)
        self.assertAlmostEqual(
            wt["total_entropy_production_rate_sigma_tot_kb_per_s"],
            wt["housekeeping_entropy_rate_sigma_hk_kb_per_s"] + wt["excess_entropy_rate_sigma_ex_kb_per_s"],
            places=4,
        )
        self.assertGreaterEqual(wt["housekeeping_fraction"], 0.50)
        self.assertGreaterEqual(wt["allosteric_stability_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores active steady state matching WT
        self.assertTrue(rescued["decomposition_equality_satisfied"])
        self.assertEqual(rescued["gating_regime"], "ACTIVE_STEADY_STATE_DISSIPATION")
        self.assertGreaterEqual(rescued["total_entropy_production_rate_sigma_tot_kb_per_s"], 2.5)
        self.assertGreaterEqual(rescued["housekeeping_entropy_rate_sigma_hk_kb_per_s"], 1.5)
        self.assertGreaterEqual(rescued["excess_entropy_rate_sigma_ex_kb_per_s"], 0.8)
        self.assertAlmostEqual(
            rescued["total_entropy_production_rate_sigma_tot_kb_per_s"],
            rescued["housekeeping_entropy_rate_sigma_hk_kb_per_s"] + rescued["excess_entropy_rate_sigma_ex_kb_per_s"],
            places=4,
        )
        self.assertGreaterEqual(rescued["housekeeping_fraction"], 0.50)
        self.assertGreaterEqual(rescued["allosteric_stability_fidelity"], 0.90)

        # F508del: Futile switching dissipation and collapsed housekeeping
        self.assertFalse(f508del["decomposition_equality_satisfied"])
        self.assertEqual(f508del["gating_regime"], "FUTILE_SWITCHING_DISSIPATION")
        self.assertLess(f508del["total_entropy_production_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(f508del["housekeeping_entropy_rate_sigma_hk_kb_per_s"], 0.1)
        self.assertLess(f508del["housekeeping_fraction"], 0.25)
        self.assertLess(f508del["allosteric_stability_fidelity"], 0.15)

        # G551D: Futile switching dissipation and collapsed housekeeping
        self.assertFalse(g551d["decomposition_equality_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_SWITCHING_DISSIPATION")
        self.assertLess(g551d["total_entropy_production_rate_sigma_tot_kb_per_s"], 0.5)
        self.assertLess(g551d["housekeeping_entropy_rate_sigma_hk_kb_per_s"], 0.1)
        self.assertLess(g551d["housekeeping_fraction"], 0.25)
        self.assertLess(g551d["allosteric_stability_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["decomposition_equality_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_rate_sigma_tot_kb_per_s"], 0.0)
        self.assertEqual(w1282x["housekeeping_entropy_rate_sigma_hk_kb_per_s"], 0.0)
        self.assertEqual(w1282x["excess_entropy_rate_sigma_ex_kb_per_s"], 0.0)
        self.assertEqual(w1282x["housekeeping_fraction"], 0.0)
        self.assertEqual(w1282x["allosteric_stability_fidelity"], 0.0)

    def test_formal_proof_dag_entropy_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        expected_benchmarks = [
            "anthropics/fermats-last-theorem",
            "anthropics/formal-math/zeta23",
            "anthropics/formal-math/3sum-apsp",
        ]
        for c_id in expected_benchmarks:
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_entropy_parity_satisfied"])
            self.assertGreaterEqual(c_prof["proof_total_entropy_rate"], 4.0)
            self.assertGreaterEqual(c_prof["proof_housekeeping_rate"], 2.5)

    def test_ah10_negative_controls(self) -> None:
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_esposito_van_den_broeck_entropy(n_gating_regimes=3)


if __name__ == "__main__":
    unittest.main()
