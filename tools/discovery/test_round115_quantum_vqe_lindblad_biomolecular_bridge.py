"""Unit tests for Round 115: Cross-Field Quantum Computing Bridge — VQE Rayleigh-Ritz & Lindblad Decoherence."""

import math
import unittest
from tools.discovery.round115_quantum_vqe_lindblad_biomolecular_bridge import (
    CITATIONS,
    HYPOTHESIS,
    PARAMETER_PROVENANCE,
    PRIMARY_DOI,
    PROFILES,
    ROUND_ID,
    SECONDARY_DOIS,
    TASK_ID,
    TITLE,
    QuantumAllostericProfile,
    compute_von_neumann_entanglement_entropy,
    evaluate_vqe_ansatz_energy,
    run_falsification,
    run_falsification_checks,
    run_round115_verification,
    simulate_quantum_vqe_lindblad,
)


class TestRound115QuantumVqeLindbladBiomolecularBridge(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round115")
        self.assertEqual(TASK_ID, "T-0750")
        self.assertEqual(PRIMARY_DOI, "10.1038/ncomms5213")
        self.assertEqual(
            SECONDARY_DOIS,
            ["10.1007/BF01608499", "10.1038/nature23879"],
        )
        self.assertIn("Variational Quantum Eigensolver", TITLE)
        self.assertIn("Lindblad", TITLE)
        self.assertIn("Rayleigh-Ritz", HYPOTHESIS)
        self.assertEqual(len(CITATIONS), 3)
        expected_dois = {
            "10.1038/ncomms5213",
            "10.1007/BF01608499",
            "10.1038/nature23879",
        }
        self.assertEqual({c["doi"] for c in CITATIONS}, expected_dois)
        self.assertTrue(all(c.get("scite_verified") is True for c in CITATIONS))

        ah10 = run_falsification()
        self.assertTrue(ah10["passed"])
        self.assertEqual(ah10["rule"], "AH-10")
        self.assertEqual(ah10["round_id"], ROUND_ID)
        self.assertEqual(ah10["task_id"], TASK_ID)

    def test_02_parameter_provenance_bounds(self) -> None:
        for p_name, p_meta in PARAMETER_PROVENANCE.items():
            self.assertIn("unit", p_meta, msg=f"Missing unit for {p_name}")
            self.assertIn("kind", p_meta, msg=f"Missing kind for {p_name}")
            self.assertIn("source", p_meta, msg=f"Missing source for {p_name}")
            self.assertIn("bounds", p_meta, msg=f"Missing bounds for {p_name}")
            low, high = p_meta["bounds"]
            self.assertLess(low, high)
            if "value" in p_meta:
                val = p_meta["value"]
                self.assertGreaterEqual(val, low)
                self.assertLessEqual(val, high)

    def test_03_five_cftr_quantum_profiles_inventory(self) -> None:
        expected_names = [
            "wt_cftr_coherent_vqe",
            "f508del_eti_error_mitigated",
            "g551d_ivacaftor_phase_locked",
            "f508del_uncorrected_decohered",
            "w1282x_nonsense_dephased",
        ]
        self.assertEqual([p.name for p in PROFILES], expected_names)

    def test_04_all_falsification_checks_and_verification_contract(self) -> None:
        report = run_falsification_checks()
        self.assertTrue(report["passed"])
        for check_name, check_val in report["checks"].items():
            self.assertTrue(check_val, msg=f"Check failed: {check_name}")

        ver = run_round115_verification()
        self.assertEqual(ver["round"], 115)
        self.assertEqual(ver["status"], "PASS")
        self.assertEqual(ver["doi"], "10.1038/ncomms5213")
        self.assertEqual(len(ver["secondary_dois"]), 2)
        self.assertEqual(len(ver["genotypes"]), 5)
        self.assertEqual(len(ver["negative_controls"]), 3)
        self.assertTrue(all(nc["rejected"] for nc in ver["negative_controls"]))
        self.assertTrue(ver["summary_metrics"]["parametric_sensitivity_25pct_passed"])

    def test_05_vqe_rayleigh_ritz_and_entanglement_entropy_invariants(self) -> None:
        for p in PROFILES:
            res = simulate_quantum_vqe_lindblad(p)
            self.assertAlmostEqual(
                res.vqe_optimal_energy_kb_t,
                res.exact_ground_energy_e0_kb_t,
                places=12,
            )
            self.assertGreaterEqual(res.vqe_grid_min_energy_kb_t, res.exact_ground_energy_e0_kb_t - 1e-12)
            self.assertGreater(res.von_neumann_entanglement_entropy_nats, 0.0)
            self.assertLessEqual(res.von_neumann_entanglement_entropy_nats, math.log(2.0) + 1e-12)
            self.assertGreater(res.concurrence, 0.0)
            self.assertLessEqual(res.concurrence, 1.0)
            self.assertTrue(res.is_rayleigh_ritz_satisfied)
            self.assertTrue(res.is_entanglement_entropy_non_negative)

    def test_06_lindblad_cptp_trace_purity_and_error_mitigation_ordering(self) -> None:
        results = {p.name: simulate_quantum_vqe_lindblad(p) for p in PROFILES}
        for res in results.values():
            self.assertAlmostEqual(res.density_matrix_trace, 1.0, places=12)
            self.assertGreaterEqual(res.density_matrix_purity, 0.25)
            self.assertLessEqual(res.density_matrix_purity, 1.0)
            self.assertGreaterEqual(res.lindblad_open_system_energy_kb_t, res.exact_ground_energy_e0_kb_t)
            self.assertTrue(res.is_lindblad_cptp_valid)

        wt = results["wt_cftr_coherent_vqe"]
        eti = results["f508del_eti_error_mitigated"]
        g551d = results["g551d_ivacaftor_phase_locked"]
        uncorr = results["f508del_uncorrected_decohered"]
        w1282x = results["w1282x_nonsense_dephased"]

        self.assertGreater(wt.ground_state_fidelity, eti.ground_state_fidelity)
        self.assertGreater(eti.ground_state_fidelity, g551d.ground_state_fidelity)
        self.assertGreater(g551d.ground_state_fidelity, uncorr.ground_state_fidelity)
        self.assertGreater(uncorr.ground_state_fidelity, w1282x.ground_state_fidelity)
        self.assertGreater(eti.ground_state_fidelity, 0.92)

    def test_07_three_poisoned_negative_controls_fail_strictly(self) -> None:
        p_rr = QuantumAllostericProfile(
            name="bad_rr",
            description="Rayleigh-Ritz violation",
            delta_eff_kb_t=1.5,
            j_cpl_kb_t=2.0,
            gamma_phi_s_inv=0.2,
            tau_s=0.25,
            poison_mode="rayleigh_ritz_violation",
        )
        res_rr = simulate_quantum_vqe_lindblad(p_rr)
        self.assertFalse(res_rr.is_rayleigh_ritz_satisfied)
        self.assertLess(res_rr.rayleigh_ritz_min_slack_kb_t, -1.0)

        p_cptp = QuantumAllostericProfile(
            name="bad_cptp",
            description="Purity trace violation",
            delta_eff_kb_t=1.5,
            j_cpl_kb_t=2.0,
            gamma_phi_s_inv=0.2,
            tau_s=0.25,
            poison_mode="purity_trace_violation",
        )
        res_cptp = simulate_quantum_vqe_lindblad(p_cptp)
        self.assertFalse(res_cptp.is_lindblad_cptp_valid)
        self.assertGreater(res_cptp.density_matrix_purity, 1.0)

        p_ent = QuantumAllostericProfile(
            name="bad_ent",
            description="Negative entanglement entropy violation",
            delta_eff_kb_t=1.5,
            j_cpl_kb_t=2.0,
            gamma_phi_s_inv=0.2,
            tau_s=0.25,
            poison_mode="negative_entanglement_entropy",
        )
        res_ent = simulate_quantum_vqe_lindblad(p_ent)
        self.assertFalse(res_ent.is_entanglement_entropy_non_negative)
        self.assertLess(res_ent.von_neumann_entanglement_entropy_nats, 0.0)
        self.assertAlmostEqual(compute_von_neumann_entanglement_entropy(0.0), 0.0, places=12)
        self.assertAlmostEqual(
            compute_von_neumann_entanglement_entropy(0.5 * math.pi),
            math.log(2.0),
            places=12,
        )
        self.assertLess(evaluate_vqe_ansatz_energy(0.5, 2.0, 3.0), 0.0)


if __name__ == "__main__":
    unittest.main()
