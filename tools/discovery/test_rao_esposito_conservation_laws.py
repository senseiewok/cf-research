"""Unit tests for Round 62 (T-0629): Rao-Esposito Conservation Laws and Topological Decomposition of Dissipation in CFTR Channel Networks in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevX.6.041064 (Rao R, Esposito M: Nonequilibrium thermodynamics of chemical reaction networks: Wisdom from stochastic thermodynamics, Phys. Rev. X 6, 041064, 2016)
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
    evaluate_rao_esposito_conservation_laws,
)


class RaoEspositoConservationLawsTests(unittest.TestCase):
    """Test suite verifying Round 62 Rao-Esposito conservation laws and topological moiety decomposition for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_rao_esposito_conservation_laws()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0629")
        self.assertTrue(self.report["rao_esposito_conservation_bridge_verified"])
        self.assertTrue(self.report["conservation_laws_verified"])
        self.assertTrue(self.report["moiety_decomposition_verified"])
        self.assertTrue(self.report["emergent_force_stratified"])
        self.assertTrue(self.report["proof_dag_conservation_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_stoichiometric_conservation_and_moiety_decomposition(self) -> None:
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

        # WT CFTR: Strict stoichiometric coupling, unbroken moiety conservation >= 0.95, high emergent force >= 18 kJ/mol
        self.assertTrue(wt["stoichiometric_conservation_satisfied"])
        self.assertEqual(wt["gating_regime"], "STRICT_STOICHIOMETRIC_COUPLING")
        self.assertGreaterEqual(wt["emergent_thermodynamic_force_kj_mol"], 18.0)
        self.assertGreater(wt["internal_dissipation_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(wt["moiety_conservation_ratio"], 0.95)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(wt["broken_moiety_leak_rate_kcal_mol_s"], 0.1)

        # Rescued CFTR (Ivacaftor / Trikafta): Restores stoichiometric coupling and moiety conservation >= 0.90
        self.assertTrue(rescued["stoichiometric_conservation_satisfied"])
        self.assertEqual(rescued["gating_regime"], "STRICT_STOICHIOMETRIC_COUPLING")
        self.assertGreaterEqual(rescued["emergent_thermodynamic_force_kj_mol"], 18.0)
        self.assertGreater(rescued["internal_dissipation_rate_kcal_mol_s"], 3.0)
        self.assertGreaterEqual(rescued["moiety_conservation_ratio"], 0.90)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertLessEqual(rescued["broken_moiety_leak_rate_kcal_mol_s"], 0.2)

        # F508del: Broken moiety leak, uncoupled dissipation
        self.assertFalse(f508del["stoichiometric_conservation_satisfied"])
        self.assertEqual(f508del["gating_regime"], "BROKEN_MOIETY_LEAK")
        self.assertLess(f508del["emergent_thermodynamic_force_kj_mol"], 5.0)
        self.assertLess(f508del["moiety_conservation_ratio"], 0.25)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(f508del["broken_moiety_leak_rate_kcal_mol_s"], 0.5)

        # G551D: Catalytic decoupler, broken stoichiometric coupling
        self.assertFalse(g551d["stoichiometric_conservation_satisfied"])
        self.assertEqual(g551d["gating_regime"], "BROKEN_MOIETY_LEAK")
        self.assertLess(g551d["emergent_thermodynamic_force_kj_mol"], 5.0)
        self.assertLess(g551d["moiety_conservation_ratio"], 0.20)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)
        self.assertGreater(g551d["broken_moiety_leak_rate_kcal_mol_s"], 0.5)

        # W1282X: Passive equilibrium baseline
        self.assertFalse(w1282x["stoichiometric_conservation_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["emergent_thermodynamic_force_kj_mol"], 0.0)
        self.assertEqual(w1282x["internal_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["broken_moiety_leak_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_formal_proof_dag_conservation_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        for c_id in ("anthropics/fermats-last-theorem", "anthropics/formal-math/zeta23", "anthropics/formal-math/3sum-apsp"):
            self.assertIn(c_id, challenges)
            c_prof = challenges[c_id]
            self.assertTrue(c_prof["proof_conservation_parity_satisfied"])
            self.assertGreater(c_prof["proof_emergent_conservation_force"], 10.0)
            self.assertGreaterEqual(c_prof["proof_moiety_conservation_ratio"], 0.85)

    def test_negative_controls(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_rao_esposito_conservation_laws(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
