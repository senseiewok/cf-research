"""Unit tests for Round 77 (T-0665): Van Vu-Hasegawa Geometrical Bounds of Irreversibility and Dynamical Activity in CFTR Catalytic Assemblies in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.126.010601 (Van Vu T, Hasegawa Y: Geometrical bounds of the irreversibility in Markovian systems, Phys. Rev. Lett. 126, 010601, 2021)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis, J. Gen. Physiol. 2006)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp, formal_math_eth
"""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_van_vu_hasegawa_geometrical_irreversibility,
)


class VanVuHasegawaGeometricalIrreversibilityTests(unittest.TestCase):
    """Test suite verifying Round 77 Van Vu-Hasegawa geometrical irreversibility & dynamical activity bounds."""

    def setUp(self) -> None:
        self.report = evaluate_van_vu_hasegawa_geometrical_irreversibility()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0665")
        self.assertTrue(self.report["van_vu_hasegawa_geometrical_irreversibility_bounds_bridge_verified"])
        self.assertTrue(self.report["geometrical_irreversibility_bound_verified"])
        self.assertTrue(self.report["dynamical_activity_scaling_verified"])
        self.assertTrue(self.report["trikafta_geometric_restoration_verified"])
        self.assertTrue(self.report["proof_dag_geometrical_irreversibility_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_gating_profiles_and_irreversibility_bounds(self) -> None:
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

        # WT CFTR: Van Vu-Hasegawa geometric bound:
        # Sigma >= 2 * L^2 / K and eta_geom = 2 * L^2 / (Sigma * K) <= 1.0
        self.assertTrue(wt["bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "GEOMETRICAL_IRREVERSIBILITY_OPTIMAL")
        self.assertGreaterEqual(wt["total_entropy_production_sigma_kb"], wt["van_vu_hasegawa_bound"])
        self.assertLessEqual(wt["geometrical_efficiency_ratio"], 1.0)
        self.assertGreater(wt["geometrical_efficiency_ratio"], 0.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertGreaterEqual(wt["total_jump_activity_k"], 15.0)
        self.assertGreaterEqual(wt["total_variation_distance_l"], 0.60)

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores geometric irreversibility bound
        self.assertTrue(rescued["bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "GEOMETRICAL_IRREVERSIBILITY_OPTIMAL")
        self.assertGreaterEqual(rescued["total_entropy_production_sigma_kb"], rescued["van_vu_hasegawa_bound"])
        self.assertLessEqual(rescued["geometrical_efficiency_ratio"], 1.0)
        self.assertGreater(rescued["geometrical_efficiency_ratio"], 0.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.85)

        # Uncorrected F508del (PDB 5UAK): Kinetic arrest & geometric activity collapse
        self.assertFalse(f508del["bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "KINETIC_ARREST_GEOMETRIC_COLLAPSE")
        self.assertLess(f508del["total_jump_activity_k"], 1.0)
        self.assertLess(f508del["total_variation_distance_l"], 0.10)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.20)

        # Uncoupled G551D: Futile uncoupled dissipation without pore displacement
        self.assertFalse(g551d["bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_UNCOUPLED_GEOMETRIC_COLLAPSE")
        self.assertLess(g551d["total_variation_distance_l"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)

        # W1282X nonsense baseline: Passive equilibrium apo
        self.assertFalse(w1282x["bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_entropy_production_sigma_kb"], 0.0)
        self.assertEqual(w1282x["total_variation_distance_l"], 0.0)
        self.assertEqual(w1282x["total_jump_activity_k"], 0.0)

    def test_challenge_profiles_formal_proof_dag_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertEqual(len(challenges), 4)
        for c_id in ["fermats_last_theorem", "formal_math_zeta23", "formal_math_3sum_apsp", "formal_math_eth"]:
            self.assertIn(c_id, challenges)
            c = challenges[c_id]
            self.assertTrue(c["parity_verified"])
            self.assertTrue(c["bound_satisfied"])
            self.assertGreater(c["sigma_proof_entropy"], 0.0)
            self.assertGreater(c["total_jump_activity_k"], 0.0)
            self.assertGreater(c["total_variation_distance_l"], 0.0)
            self.assertGreaterEqual(c["sigma_proof_entropy"], c["van_vu_hasegawa_bound"])
            self.assertLessEqual(c["geometrical_efficiency_ratio"], 1.0)

    def test_n_gating_regimes_validation(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_van_vu_hasegawa_geometrical_irreversibility(n_gating_regimes=4)
        with self.assertRaises(ValueError):
            evaluate_van_vu_hasegawa_geometrical_irreversibility(n_gating_regimes=6)


if __name__ == "__main__":
    unittest.main()
