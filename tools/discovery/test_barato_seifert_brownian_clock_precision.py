"""Unit tests for Round 82 (T-0675): Barato-Seifert Cost of Precision & Biomolecular Clock Efficiency in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevX.6.041053 (Barato AC, Seifert U: Cost and Precision of Brownian Clocks, Phys. Rev. X 6, 041053, 2016)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis, J. Gen. Physiol. 2006)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM, 5UAK, 8EIQ
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_barato_seifert_brownian_clock_precision,
)


class BaratoSeifertBrownianClockPrecisionTests(unittest.TestCase):
    """Test suite verifying Round 82 Barato-Seifert Brownian clock cost-precision bounds."""

    def setUp(self) -> None:
        self.report = evaluate_barato_seifert_brownian_clock_precision()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0675")
        self.assertTrue(self.report["barato_seifert_brownian_clock_precision_bridge_verified"])
        self.assertTrue(self.report["clock_cost_precision_bound_verified"])
        self.assertTrue(self.report["optimal_brownian_clock_verified"])
        self.assertTrue(self.report["trikafta_clock_restoration_verified"])
        self.assertTrue(self.report["proof_dag_clock_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_clock_bounds(self) -> None:
        profiles = self.report["gating_profiles"]
        for gid in ["wt_cftr", "rescued_cftr_ivacaftor", "f508del_uncorrected", "g551d_uncoupled_cftr", "w1282x_nonsense_baseline"]:
            self.assertIn(gid, profiles)

        wt = profiles["wt_cftr"]
        rescued = profiles["rescued_cftr_ivacaftor"]
        f508del = profiles["f508del_uncorrected"]
        g551d = profiles["g551d_uncoupled_cftr"]
        w1282x = profiles["w1282x_nonsense_baseline"]

        self.assertEqual(wt["gating_regime"], "AUTONOMOUS_BROWNIAN_CLOCK_OPTIMAL")
        self.assertTrue(wt["bound_satisfied"])
        self.assertLessEqual(wt["observed_clock_precision_n"], wt["theoretical_max_precision_n_max"])
        self.assertGreaterEqual(wt["clock_efficiency_eta"], 0.85)

        self.assertEqual(rescued["gating_regime"], "RESCUED_BROWNIAN_CLOCK_RESTORED")
        self.assertTrue(rescued["bound_satisfied"])
        self.assertLessEqual(rescued["observed_clock_precision_n"], rescued["theoretical_max_precision_n_max"])
        self.assertGreaterEqual(rescued["clock_efficiency_eta"], 0.85)

        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_CLOCK_DECOHERENCE")
        self.assertFalse(f508del["bound_satisfied"])
        self.assertLess(f508del["clock_efficiency_eta"], 0.25)

        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_CLOCK_COLLAPSE")
        self.assertFalse(g551d["bound_satisfied"])
        self.assertLess(g551d["clock_efficiency_eta"], 0.20)

        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertFalse(w1282x["bound_satisfied"])
        self.assertEqual(w1282x["clock_efficiency_eta"], 0.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["clock_bound_satisfied"])
            self.assertGreaterEqual(c["proof_clock_efficiency_eta"], 0.85)

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_barato_seifert_brownian_clock_precision(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
