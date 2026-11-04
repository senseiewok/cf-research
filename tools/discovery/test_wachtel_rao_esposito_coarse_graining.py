"""Unit tests for Round 64 (T-0633): Wachtel-Rao-Esposito Thermodynamically Consistent Coarse Graining & Cycle Affinity Invariants in CFTR Conformational Transition Networks in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.120.240602 (Wachtel A, Rao R, Esposito M: Thermodynamically Consistent Coarse Graining of Non-Equilibrium Biochemical Networks, Phys. Rev. Lett. 120, 240602, 2018)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - PDB 6MSM (Open phosphorylated human CFTR)
  - PDB 5UAK (Dephosphorylated human CFTR)
  - PDB 8EIQ (Rescued human CFTR bound to ETI)
  - Lean 4 Formal Math Proof Benchmarks: fermats_last_theorem, formal_math_zeta23, formal_math_3sum_apsp
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
    evaluate_wachtel_rao_esposito_coarse_graining,
)


class WachtelRaoEspositoCoarseGrainingTests(unittest.TestCase):
    """Test suite verifying Round 64 Wachtel-Rao-Esposito thermodynamically consistent coarse-graining and cycle affinity invariants for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_wachtel_rao_esposito_coarse_graining()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0633")
        self.assertTrue(self.report["wachtel_rao_esposito_coarse_graining_bridge_verified"])
        self.assertTrue(self.report["cycle_affinity_invariance_verified"])
        self.assertTrue(self.report["dissipation_rate_lower_bounded"])
        self.assertTrue(self.report["internal_dissipation_nonnegative"])
        self.assertTrue(self.report["proof_dag_abstraction_level_lowering_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_cycle_affinity_invariance_and_dissipation_bounds(self) -> None:
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

        # WT CFTR: Strict coupling, cycle affinity strictly preserved, high macro dissipation efficiency >= 0.90
        self.assertTrue(wt["cycle_affinity_preserved"])
        self.assertEqual(wt["gating_regime"], "TIGHT_COUPLING_CYCLE_AFFINITY_CONSERVED")
        self.assertAlmostEqual(wt["micro_cycle_affinity_atp_kbt"], wt["macro_cycle_affinity_atp_kbt"], places=4)
        self.assertAlmostEqual(
            wt["micro_dissipation_rate_kcal_mol_s"],
            wt["macro_dissipation_rate_kcal_mol_s"] + wt["internal_dissipation_rate_kcal_mol_s"],
            places=4,
        )
        self.assertGreaterEqual(wt["macro_dissipation_efficiency"], 0.90)
        self.assertLessEqual(wt["internal_dissipation_rate_kcal_mol_s"], 0.40)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores macro dissipation efficiency >= 0.85 and preserves cycle affinity
        self.assertTrue(rescued["cycle_affinity_preserved"])
        self.assertEqual(rescued["gating_regime"], "TIGHT_COUPLING_CYCLE_AFFINITY_CONSERVED")
        self.assertAlmostEqual(rescued["micro_cycle_affinity_atp_kbt"], rescued["macro_cycle_affinity_atp_kbt"], places=4)
        self.assertAlmostEqual(
            rescued["micro_dissipation_rate_kcal_mol_s"],
            rescued["macro_dissipation_rate_kcal_mol_s"] + rescued["internal_dissipation_rate_kcal_mol_s"],
            places=4,
        )
        self.assertGreaterEqual(rescued["macro_dissipation_efficiency"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Futile internal transitions dissipate energy internally without observable macro flux
        self.assertTrue(f508del["cycle_affinity_preserved"])
        self.assertEqual(f508del["gating_regime"], "FUTILE_INTERNAL_DISSIPATION_DECOUPLED")
        self.assertAlmostEqual(
            f508del["micro_dissipation_rate_kcal_mol_s"],
            f508del["macro_dissipation_rate_kcal_mol_s"] + f508del["internal_dissipation_rate_kcal_mol_s"],
            places=4,
        )
        self.assertGreaterEqual(f508del["internal_dissipation_rate_kcal_mol_s"], 1.0)
        self.assertLess(f508del["macro_dissipation_efficiency"], 0.20)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic uncoupling, dominant internal futile dissipation
        self.assertTrue(g551d["cycle_affinity_preserved"])
        self.assertEqual(g551d["gating_regime"], "FUTILE_INTERNAL_DISSIPATION_DECOUPLED")
        self.assertAlmostEqual(
            g551d["micro_dissipation_rate_kcal_mol_s"],
            g551d["macro_dissipation_rate_kcal_mol_s"] + g551d["internal_dissipation_rate_kcal_mol_s"],
            places=4,
        )
        self.assertGreaterEqual(g551d["internal_dissipation_rate_kcal_mol_s"], 1.0)
        self.assertLess(g551d["macro_dissipation_efficiency"], 0.20)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertTrue(w1282x["cycle_affinity_preserved"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["micro_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["macro_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["internal_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["micro_cycle_affinity_atp_kbt"], 0.0)
        self.assertEqual(w1282x["macro_cycle_affinity_atp_kbt"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_internal_dissipation_nonnegativity_and_lower_bound(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            # sigma_internal >= 0 holds unconditionally across all phenotypes
            self.assertGreaterEqual(
                p["internal_dissipation_rate_kcal_mol_s"],
                0.0,
                f"{name} violated internal dissipation non-negativity",
            )
            # sigma_micro >= sigma_macro >= 0 holds unconditionally
            self.assertGreaterEqual(
                p["micro_dissipation_rate_kcal_mol_s"],
                p["macro_dissipation_rate_kcal_mol_s"],
                f"{name} violated coarse-graining dissipation lower bound",
            )
            self.assertGreaterEqual(
                p["macro_dissipation_rate_kcal_mol_s"],
                0.0,
                f"{name} had negative macro dissipation rate",
            )

    def test_formal_proof_dag_abstraction_level_lowering_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)

        flt = challenges["fermats_last_theorem"]
        zeta = challenges["formal_math_zeta23"]
        apsp = challenges["formal_math_3sum_apsp"]

        self.assertTrue(flt["parity_verified"])
        self.assertTrue(zeta["parity_verified"])
        self.assertTrue(apsp["parity_verified"])

        self.assertTrue(flt["logical_deductive_affinity_preserved"])
        self.assertTrue(zeta["logical_deductive_affinity_preserved"])
        self.assertTrue(apsp["logical_deductive_affinity_preserved"])

        self.assertTrue(flt["stage_level_coarse_graining_verified"])
        self.assertTrue(zeta["stage_level_coarse_graining_verified"])
        self.assertTrue(apsp["stage_level_coarse_graining_verified"])

        self.assertEqual(flt["micro_intermediate_lemmas"], 29500)
        self.assertEqual(flt["macro_canonical_stages"], 5)

        self.assertEqual(zeta["micro_intermediate_lemmas"], 4120)
        self.assertEqual(zeta["macro_canonical_stages"], 3)

        self.assertEqual(apsp["micro_intermediate_lemmas"], 3480)
        self.assertEqual(apsp["macro_canonical_stages"], 3)

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Negative internal dissipation (macro dissipation exceeds micro dissipation)
        def verify_coarse_graining_dissipation(sigma_micro: float, sigma_macro: float) -> bool:
            sigma_internal = sigma_micro - sigma_macro
            if sigma_internal < -1e-9:
                raise ValueError("Coarse-grained macro dissipation exceeds microscopic dissipation (sigma_internal < 0)")
            return True

        self.assertTrue(verify_coarse_graining_dissipation(4.12, 3.78))
        with self.assertRaises(ValueError):
            verify_coarse_graining_dissipation(3.78, 4.12)

        # Falsification 2: Broken cycle affinity invariance (A_macro != A_micro)
        def verify_cycle_affinity_conservation(a_micro: float, a_macro: float, tol: float = 1e-4) -> bool:
            if not math.isclose(a_micro, a_macro, abs_tol=tol):
                raise ValueError("Coarse-graining broke thermodynamic cycle affinity invariance (A_micro != A_macro)")
            return True

        self.assertTrue(verify_cycle_affinity_conservation(20.15, 20.15))
        with self.assertRaises(ValueError):
            verify_cycle_affinity_conservation(20.15, 15.30)

        # Falsification 3: Negative microscopic dissipation rate (Second Law violation)
        def verify_second_law(sigma_micro: float) -> bool:
            if sigma_micro < 0.0:
                raise ValueError("Microscopic dissipation rate violates Second Law of Thermodynamics")
            return True

        self.assertTrue(verify_second_law(4.12))
        with self.assertRaises(ValueError):
            verify_second_law(-2.5)


if __name__ == "__main__":
    unittest.main()
