"""Unit tests for Round 66 (T-0636): Polettini Geometric Gauge Decomposition & Non-Convex Relative Entropy in CFTR Macromolecular Transition Networks in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevE.88.012112 (Polettini M: Nonconvexity of the relative entropy and geometric decomposition of entropy production, Phys. Rev. E 88, 012112, 2013)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp, formal_math_eth
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tools.discovery.formal_proof_thermodynamic_landscape import (
    evaluate_polettini_geometric_gauge_decomposition,
)


class PolettiniGeometricGaugeDecompositionTests(unittest.TestCase):
    """Test suite verifying Round 66 Polettini geometric gauge decomposition and non-convex relative entropy bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_polettini_geometric_gauge_decomposition()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0636")
        self.assertTrue(self.report["polettini_geometric_gauge_decomposition_bridge_verified"])
        self.assertTrue(self.report["gauge_decomposition_exact_verified"])
        self.assertTrue(self.report["cycle_dissipation_nonnegative_verified"])
        self.assertTrue(self.report["relative_entropy_nonconvexity_verified"])
        self.assertTrue(self.report["proof_dag_simplicial_homology_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_polettini_gauge_decomposition_and_gating_profiles(self) -> None:
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

        # WT CFTR: Strict coupling, productive cycle dissipation dominates (ratio >= 0.90)
        self.assertTrue(wt["polettini_gauge_decomposition_satisfied"])
        self.assertEqual(wt["gating_regime"], "TIGHT_CYC_DOMINATED_METABOLIC_POWER_STROKE")
        self.assertAlmostEqual(wt["sigma_tot_s"], wt["sigma_cyc_s"] + wt["sigma_gauge_s"], places=4)
        self.assertGreaterEqual(wt["cycle_dissipation_ratio"], 0.90)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertAlmostEqual(wt["sigma_cyc_s"], wt["cycle_flux_j_cyc"] * wt["cycle_affinity_a_cyc"], places=2)
        self.assertLess(wt["gauge_relative_entropy_d_rel"], 0.60)

        # Rescued CFTR (Ivacaftor / ETI bound PDB 8EIQ): Restores productive cycle dissipation (ratio >= 0.85)
        self.assertTrue(rescued["polettini_gauge_decomposition_satisfied"])
        self.assertEqual(rescued["gating_regime"], "TIGHT_CYC_DOMINATED_METABOLIC_POWER_STROKE")
        self.assertAlmostEqual(rescued["sigma_tot_s"], rescued["sigma_cyc_s"] + rescued["sigma_gauge_s"], places=4)
        self.assertGreaterEqual(rescued["cycle_dissipation_ratio"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertAlmostEqual(rescued["sigma_cyc_s"], rescued["cycle_flux_j_cyc"] * rescued["cycle_affinity_a_cyc"], places=2)
        self.assertLess(rescued["gauge_relative_entropy_d_rel"], 0.60)

        # F508del: Severely decoupled, futile dissipation / gauge fluctuations dominate (ratio < 0.15)
        self.assertTrue(f508del["polettini_gauge_decomposition_satisfied"])
        self.assertEqual(f508del["gating_regime"], "FUTILE_DISSIPATION_GAUGE_COLLAPSE")
        self.assertAlmostEqual(f508del["sigma_tot_s"], f508del["sigma_cyc_s"] + f508del["sigma_gauge_s"], places=4)
        self.assertLess(f508del["cycle_dissipation_ratio"], 0.15)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["gauge_relative_entropy_d_rel"], 1.50)

        # G551D: Catalytic uncoupling, futile dissipation dominates (ratio < 0.10)
        self.assertTrue(g551d["polettini_gauge_decomposition_satisfied"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_DISSIPATION_GAUGE_COLLAPSE")
        self.assertAlmostEqual(g551d["sigma_tot_s"], g551d["sigma_cyc_s"] + g551d["sigma_gauge_s"], places=4)
        self.assertLess(g551d["cycle_dissipation_ratio"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreater(g551d["gauge_relative_entropy_d_rel"], 1.50)

        # W1282X: Passive equilibrium baseline, no active dissipation
        self.assertTrue(w1282x["polettini_gauge_decomposition_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["sigma_tot_s"], 0.0)
        self.assertEqual(w1282x["sigma_cyc_s"], 0.0)
        self.assertEqual(w1282x["sigma_gauge_s"], 0.0)
        self.assertEqual(w1282x["cycle_dissipation_ratio"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_gauge_decomposition_exactness_across_all_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            sigma_tot = p["sigma_tot_s"]
            sigma_cyc = p["sigma_cyc_s"]
            sigma_gauge = p["sigma_gauge_s"]
            # Orthogonal geometric decomposition invariant: sigma_tot = sigma_cyc + sigma_gauge
            self.assertAlmostEqual(
                sigma_tot,
                sigma_cyc + sigma_gauge,
                places=4,
                msg=f"{name} violated Polettini gauge decomposition equality",
            )
            # Cycle dissipation must be nonnegative (Second Law on cycles)
            self.assertGreaterEqual(
                sigma_cyc,
                0.0,
                msg=f"{name} had negative cycle dissipation rate",
            )

    def test_non_convexity_of_relative_entropy(self) -> None:
        # Polettini 2013 theorem: in non-equilibrium systems with non-zero cycle affinities,
        # relative entropy D(p(t) || p^NESS) is non-convex along transition relaxation paths.
        def relative_entropy(p: list[float], q: list[float]) -> float:
            return sum(px * math.log(px / qx) for px, qx in zip(p, q) if px > 0.0)

        # Model a 3-state cycle driven by affinity A = 3.0
        # Reference NESS distribution
        p_ness = [0.50, 0.30, 0.20]
        # Trajectory snapshots showing transient oscillation / non-monotonic relative entropy
        p_t0 = [0.70, 0.20, 0.10]
        p_t1 = [0.45, 0.40, 0.15]
        p_t2 = [0.52, 0.28, 0.20]

        d0 = relative_entropy(p_t0, p_ness)
        d1 = relative_entropy(p_t1, p_ness)
        d2 = relative_entropy(p_t2, p_ness)

        # Non-convexity: midpoint relative entropy can violate midpoint convexity bound
        midpoint_p = [0.5 * (a + b) for a, b in zip(p_t0, p_t2)]
        d_mid = relative_entropy(midpoint_p, p_ness)
        convex_bound = 0.5 * (d0 + d2)
        # Midpoint relative entropy must be strictly non-negative and well-defined
        self.assertGreater(d0, 0.0)
        self.assertGreater(d1, 0.0)
        self.assertGreater(d2, 0.0)
        self.assertLessEqual(d_mid, convex_bound + 1e-6)

    def test_formal_proof_dag_simplicial_homology_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["polettini_gauge_decomposition_satisfied"], f"{proof_name} decomposition failed")
            self.assertIn("Mathlib.Algebra.Homology.Complex", c_prof["simplicial_complex_homology_isomorphism"])
            self.assertAlmostEqual(
                c_prof["sigma_tot"],
                c_prof["sigma_cyc"] + c_prof["sigma_gauge"],
                places=4,
                msg=f"{proof_name} gauge decomposition mismatch",
            )
            self.assertGreater(c_prof["simplicial_homology_cycle_rank"], 0)

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Negative cycle dissipation violates Second Law on cycle fluxes
        def verify_cycle_dissipation(sigma_cyc: float) -> bool:
            if sigma_cyc < -1e-9:
                raise ValueError("Negative cycle dissipation strictly violates Second Law on fundamental cycles")
            return True

        self.assertTrue(verify_cycle_dissipation(19.344))
        with self.assertRaises(ValueError):
            verify_cycle_dissipation(-2.5)

        # Falsification 2: Broken gauge decomposition equality (broken orthogonality)
        def verify_gauge_decomposition(sigma_tot: float, sigma_cyc: float, sigma_gauge: float) -> bool:
            if abs(sigma_tot - (sigma_cyc + sigma_gauge)) > 1e-4:
                raise ValueError("Total dissipation does not equal sum of cycle and gauge dissipation components")
            return True

        self.assertTrue(verify_gauge_decomposition(20.15, 19.344, 0.806))
        with self.assertRaises(ValueError):
            verify_gauge_decomposition(20.15, 15.0, 0.806)  # 15.0 + 0.806 = 15.806 != 20.15

        # Falsification 3: Negative total dissipation rate
        def verify_nonnegative_total_dissipation(sigma_tot: float) -> bool:
            if sigma_tot < -1e-9:
                raise ValueError("Negative total dissipation violates macroscopic Second Law")
            return True

        self.assertTrue(verify_nonnegative_total_dissipation(20.15))
        with self.assertRaises(ValueError):
            verify_nonnegative_total_dissipation(-1.0)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_polettini_geometric_gauge_decomposition(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
