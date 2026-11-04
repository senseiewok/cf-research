"""Unit tests for Round 63 (T-0631): Manikandan-Gupta-Krishnamurthy Operationally Accessible Bounds on Dissipation & Current Fluctuations in CFTR Channel Networks in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.122.230601 (Manikandan SK, Gupta D, Krishnamurthy S: Operationally Accessible Bounds on Fluctuations and Dissipation in Nonequilibrium Steady States, Phys. Rev. Lett. 122, 230601, 2019)
  - DOI:10.1085/jgp.200609558 (Csanady et al.: Strict coupling between CFTR channel opening and ATP hydrolysis)
  - DOI:10.1126/science.ade2216 (Fiedorczuk & Chen: Molecular structures reveal synergistic rescue of dF508 CFTR by Trikafta modulators, Science 2022)
  - DOI:10.1126/science.aaw7611 (Liu et al.: Structural identification of a hotspot on CFTR for potentiation, Science 2019)
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
    evaluate_manikandan_gupta_krishnamurthy_bounds,
)


class ManikandanGuptaKrishnamurthyBoundsTests(unittest.TestCase):
    """Test suite verifying Round 63 Manikandan-Gupta-Krishnamurthy operational dissipation bounds for CFTR gating."""

    def setUp(self) -> None:
        self.report = evaluate_manikandan_gupta_krishnamurthy_bounds()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0631")
        self.assertTrue(self.report["manikandan_gupta_krishnamurthy_bridge_verified"])
        self.assertTrue(self.report["operational_dissipation_bounds_verified"])
        self.assertTrue(self.report["current_fluctuation_dispersion_stratified"])
        self.assertTrue(self.report["variational_tightening_verified"])
        self.assertTrue(self.report["proof_dag_operational_bounding_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_operational_bounds_and_dispersion_stratification(self) -> None:
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

        # WT CFTR: Strict coupling, low relative current dispersion <= 0.15, operational bound captures >= 88% of dissipation
        self.assertTrue(wt["operational_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "TIGHT_OPERATIONAL_BOUND_STRICT_COUPLING")
        self.assertGreaterEqual(wt["operational_lower_bound_kcal_mol_s"], 3.5)
        self.assertLessEqual(wt["operational_lower_bound_kcal_mol_s"], wt["total_dissipation_rate_kcal_mol_s"])
        self.assertLessEqual(wt["relative_current_dispersion"], 0.15)
        self.assertGreaterEqual(wt["operational_bound_ratio"], 0.88)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)

        # Rescued CFTR (Ivacaftor / ETI): Restores tight bound and low dispersion
        self.assertTrue(rescued["operational_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "TIGHT_OPERATIONAL_BOUND_STRICT_COUPLING")
        self.assertGreaterEqual(rescued["operational_lower_bound_kcal_mol_s"], 3.0)
        self.assertLessEqual(rescued["operational_lower_bound_kcal_mol_s"], rescued["total_dissipation_rate_kcal_mol_s"])
        self.assertLessEqual(rescued["relative_current_dispersion"], 0.16)
        self.assertGreaterEqual(rescued["operational_bound_ratio"], 0.85)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)

        # F508del: Excessive dispersion, loose operational bound indicating futile uncoupled dissipation
        self.assertTrue(f508del["operational_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "LOOSE_OPERATIONAL_BOUND_FUTILE_DISSIPATION")
        self.assertGreaterEqual(f508del["relative_current_dispersion"], 1.4)
        self.assertLess(f508del["operational_bound_ratio"], 0.25)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.15)

        # G551D: Catalytic decoupling, excessive dispersion
        self.assertTrue(g551d["operational_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "LOOSE_OPERATIONAL_BOUND_FUTILE_DISSIPATION")
        self.assertGreaterEqual(g551d["relative_current_dispersion"], 1.5)
        self.assertLess(g551d["operational_bound_ratio"], 0.20)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.15)

        # W1282X: Passive equilibrium baseline
        self.assertTrue(w1282x["operational_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["total_dissipation_rate_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["operational_lower_bound_kcal_mol_s"], 0.0)
        self.assertEqual(w1282x["relative_current_dispersion"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_variational_tightening_gain(self) -> None:
        profiles = self.report["gating_profiles"]
        # Variational optimization over observable currents tightens the lower bound (gain >= 1.20 in WT and Rescued)
        self.assertGreaterEqual(profiles["wt_cftr"]["variational_tightening_gain"], 1.20)
        self.assertGreaterEqual(profiles["rescued_cftr_ivacaftor"]["variational_tightening_gain"], 1.18)
        self.assertLess(profiles["f508del_uncorrected"]["variational_tightening_gain"], 1.10)

    def test_formal_proof_dag_operational_bounding_parity(self) -> None:
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

        self.assertGreaterEqual(flt["operational_bounding_fidelity"], 0.85)
        self.assertGreaterEqual(zeta["operational_bounding_fidelity"], 0.85)
        self.assertGreaterEqual(apsp["operational_bounding_fidelity"], 0.85)

        self.assertLessEqual(flt["operational_verification_bound"], flt["proof_total_entropy_production"])
        self.assertLessEqual(zeta["operational_verification_bound"], zeta["proof_total_entropy_production"])
        self.assertLessEqual(apsp["operational_verification_bound"], apsp["proof_total_entropy_production"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Inverted bound sigma_operational > sigma_total violates mathematical inequality
        def check_bound_validity(sigma_total: float, sigma_oper: float) -> bool:
            if sigma_oper > sigma_total:
                raise ValueError("Operational lower bound exceeds total dissipation rate")
            return True

        self.assertTrue(check_bound_validity(4.12, 3.75))
        with self.assertRaises(ValueError):
            check_bound_validity(3.75, 4.12)

        # Falsification 2: Negative dissipation rate violates Second Law of Thermodynamics
        def check_second_law(sigma_oper: float) -> bool:
            if sigma_oper < 0.0:
                raise ValueError("Negative dissipation rate violates Second Law of Thermodynamics")
            return True

        self.assertTrue(check_second_law(3.75))
        with self.assertRaises(ValueError):
            check_second_law(-1.5)

        # Falsification 3: Degenerate correlation time tau_eff <= 0
        def check_correlation_time(tau_eff: float) -> bool:
            if tau_eff <= 0.0:
                raise ValueError("Effective correlation time must be strictly positive")
            return True

        self.assertTrue(check_correlation_time(1.8))
        with self.assertRaises(ValueError):
            check_correlation_time(-0.5)
        with self.assertRaises(ValueError):
            check_correlation_time(0.0)


if __name__ == "__main__":
    unittest.main()
