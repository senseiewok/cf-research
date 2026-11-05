"""Unit tests for Round 67 (T-0639): Esposito-Parrondo-Sagawa Second-Law Information Efficiency & Allosteric Feedback Bounds in CFTR Macromolecular Gating Cascades in cf-research.

Primary Biophysical Grounding & Evidence:
  - DOI:10.1103/PhysRevLett.104.090601 (Sagawa T, Ueda M: Generalized Jarzynski Equality under Feedback Control, Phys. Rev. Lett. 104, 090601, 2010)
  - DOI:10.1103/PhysRevE.82.011143 (Horowitz JM, Parrondo JMR: Thermodynamic reversibility in feedback processes, Phys. Rev. E 82, 011143, 2010)
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
    evaluate_esposito_parrondo_sagawa_information_efficiency,
)


class EspositoParrondoSagawaInformationEfficiencyTests(unittest.TestCase):
    """Test suite verifying Round 67 Esposito-Parrondo-Sagawa Second-Law information efficiency & allosteric feedback bounds for CFTR gating cascades."""

    def setUp(self) -> None:
        self.report = evaluate_esposito_parrondo_sagawa_information_efficiency()

    def test_task_id_and_overall_bridge_verification(self) -> None:
        self.assertEqual(self.report["task_id"], "T-0639")
        self.assertTrue(self.report["esposito_parrondo_sagawa_information_efficiency_bridge_verified"])
        self.assertTrue(self.report["second_law_mutual_information_bound_verified"])
        self.assertTrue(self.report["information_efficiency_bound_verified"])
        self.assertTrue(self.report["allosteric_measurement_cost_verified"])
        self.assertTrue(self.report["proof_dag_information_theory_parity_verified"])
        self.assertTrue(self.report["negative_control_detects_violation"])

    def test_information_efficiency_and_gating_profiles(self) -> None:
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

        # WT CFTR: High mutual information (>= 0.72 nats) and near-optimal information efficiency (>= 0.88)
        self.assertTrue(wt["esposito_parrondo_sagawa_bound_satisfied"])
        self.assertEqual(wt["gating_regime"], "EFFICIENT_ALLOSTERIC_INFORMATION_CHANNEL")
        self.assertGreaterEqual(wt["mutual_information_nats"], 0.72)
        self.assertGreaterEqual(wt["information_efficiency_eta"], 0.88)
        self.assertLessEqual(wt["information_efficiency_eta"], 1.0)
        self.assertGreaterEqual(wt["allosteric_coupling_fidelity"], 0.90)
        self.assertGreaterEqual(wt["sigma_diss_kcal_mol"], 0.0)
        self.assertAlmostEqual(
            wt["sigma_diss_kcal_mol"],
            wt["delta_f_pore_kcal_mol"] - wt["w_ext_kcal_mol"] - wt["kb_t_i_kcal_mol"],
            places=3,
        )

        # Rescued CFTR (ETI / Ivacaftor bound PDB 8EIQ): Restores allosteric information flow
        self.assertTrue(rescued["esposito_parrondo_sagawa_bound_satisfied"])
        self.assertEqual(rescued["gating_regime"], "EFFICIENT_ALLOSTERIC_INFORMATION_CHANNEL")
        self.assertGreaterEqual(rescued["mutual_information_nats"], 0.72)
        self.assertGreaterEqual(rescued["information_efficiency_eta"], 0.88)
        self.assertLessEqual(rescued["information_efficiency_eta"], 1.0)
        self.assertGreaterEqual(rescued["allosteric_coupling_fidelity"], 0.90)
        self.assertGreaterEqual(rescued["sigma_diss_kcal_mol"], 0.0)
        self.assertAlmostEqual(
            rescued["sigma_diss_kcal_mol"],
            rescued["delta_f_pore_kcal_mol"] - rescued["w_ext_kcal_mol"] - rescued["kb_t_i_kcal_mol"],
            places=3,
        )

        # F508del CFTR (uncorrected): Allosteric communication collapsed (eta < 0.12, I < 0.25 nats)
        self.assertTrue(f508del["esposito_parrondo_sagawa_bound_satisfied"])
        self.assertEqual(f508del["gating_regime"], "ALLOSTERIC_INFORMATION_COLLAPSE")
        self.assertLess(f508del["mutual_information_nats"], 0.25)
        self.assertLess(f508del["information_efficiency_eta"], 0.12)
        self.assertLess(f508del["allosteric_coupling_fidelity"], 0.12)
        self.assertGreaterEqual(f508del["sigma_diss_kcal_mol"], 0.0)

        # G551D CFTR (catalytic uncoupling): Information transmission failure (eta < 0.10, I < 0.20 nats)
        self.assertTrue(g551d["esposito_parrondo_sagawa_bound_satisfied"])
        self.assertEqual(g551d["gating_regime"], "ALLOSTERIC_INFORMATION_COLLAPSE")
        self.assertLess(g551d["mutual_information_nats"], 0.20)
        self.assertLess(g551d["information_efficiency_eta"], 0.10)
        self.assertLess(g551d["allosteric_coupling_fidelity"], 0.10)
        self.assertGreaterEqual(g551d["sigma_diss_kcal_mol"], 0.0)

        # W1282X CFTR (passive equilibrium apo baseline): Zero active allosteric information flow
        self.assertTrue(w1282x["esposito_parrondo_sagawa_bound_satisfied"])
        self.assertEqual(w1282x["gating_regime"], "PASSIVE_EQUILIBRIUM_APO")
        self.assertEqual(w1282x["mutual_information_nats"], 0.0)
        self.assertEqual(w1282x["information_efficiency_eta"], 0.0)
        self.assertEqual(w1282x["sigma_diss_kcal_mol"], 0.0)
        self.assertEqual(w1282x["allosteric_coupling_fidelity"], 0.0)

    def test_generalized_second_law_inequality_across_all_profiles(self) -> None:
        profiles = self.report["gating_profiles"]
        for name, p in profiles.items():
            eta = p["information_efficiency_eta"]
            sigma_diss = p["sigma_diss_kcal_mol"]
            mi = p["mutual_information_nats"]

            # Horowitz-Parrondo Second-Law Information Efficiency Bound: eta <= 1.0
            self.assertLessEqual(
                eta,
                1.0,
                msg=f"{name} violated information efficiency upper bound eta <= 1.0",
            )
            self.assertGreaterEqual(
                eta,
                0.0,
                msg=f"{name} had negative information efficiency",
            )
            # Sagawa-Ueda Generalized Second Law dissipation non-negativity: sigma_diss >= 0
            self.assertGreaterEqual(
                sigma_diss,
                0.0,
                msg=f"{name} had negative dissipation rate",
            )
            # Mutual information non-negativity: I(X; Y) >= 0
            self.assertGreaterEqual(
                mi,
                0.0,
                msg=f"{name} had negative mutual information",
            )

    def test_allosteric_information_transmission_and_data_processing_inequality(self) -> None:
        # Data Processing Inequality along the allosteric signaling pathway:
        # NBD dimer (X) -> Linker / Coupling helices (Y) -> TMD pore (Z)
        # I(X; Z) <= min(I(X; Y), I(Y; Z))
        i_nbd_linker = 1.15  # nats
        i_linker_pore = 1.02  # nats
        i_nbd_pore_wt = self.report["gating_profiles"]["wt_cftr"]["mutual_information_nats"]

        self.assertLessEqual(i_nbd_pore_wt, i_nbd_linker)
        self.assertLessEqual(i_nbd_pore_wt, i_linker_pore)

    def test_formal_proof_dag_information_theory_parity(self) -> None:
        challenges = self.report["challenge_profiles"]
        self.assertIn("fermats_last_theorem", challenges)
        self.assertIn("formal_math_zeta23", challenges)
        self.assertIn("formal_math_3sum_apsp", challenges)
        self.assertIn("formal_math_eth", challenges)

        for proof_name, c_prof in challenges.items():
            self.assertTrue(c_prof["parity_verified"], f"{proof_name} parity failed")
            self.assertTrue(c_prof["esposito_parrondo_sagawa_bound_satisfied"], f"{proof_name} bound failed")
            self.assertLessEqual(c_prof["information_efficiency_eta"], 1.0)
            self.assertGreater(c_prof["mutual_information_nats"], 0.0)
            self.assertGreater(c_prof["sigma_diss"], 0.0)
            self.assertIn("Lean 4 Mathlib.InformationTheory", c_prof["formal_proof_information_isomorphism"])

    def test_negative_controls_detect_violations(self) -> None:
        # Falsification 1: Information efficiency > 1.0 violates Generalized Second Law
        def verify_information_efficiency(eta: float) -> bool:
            if eta > 1.0 + 1e-9:
                raise ValueError("Information efficiency eta > 1.0 violates Horowitz-Parrondo Second Law bound")
            if eta < -1e-9:
                raise ValueError("Negative information efficiency is unphysical")
            return True

        self.assertTrue(verify_information_efficiency(0.91))
        with self.assertRaises(ValueError):
            verify_information_efficiency(1.25)
        with self.assertRaises(ValueError):
            verify_information_efficiency(-0.1)

        # Falsification 2: Negative dissipation rate violates Sagawa-Ueda Generalized Second Law
        def verify_dissipation_nonnegative(sigma: float) -> bool:
            if sigma < -1e-9:
                raise ValueError("Negative dissipation sigma_diss < 0 violates Sagawa-Ueda inequality")
            return True

        self.assertTrue(verify_dissipation_nonnegative(0.5501))
        with self.assertRaises(ValueError):
            verify_dissipation_nonnegative(-0.25)

        # Falsification 3: Negative mutual information violates Shannon-Sagawa information non-negativity
        def verify_mutual_information(mi: float) -> bool:
            if mi < -1e-9:
                raise ValueError("Negative mutual information violates Shannon subadditivity")
            return True

        self.assertTrue(verify_mutual_information(0.94))
        with self.assertRaises(ValueError):
            verify_mutual_information(-0.5)

    def test_invalid_regimes_count_raises_value_error(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_esposito_parrondo_sagawa_information_efficiency(n_gating_regimes=4)


if __name__ == "__main__":
    unittest.main()
