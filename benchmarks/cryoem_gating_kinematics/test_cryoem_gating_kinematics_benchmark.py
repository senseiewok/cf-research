#!/usr/bin/env python3
"""Unit tests for the Modulator-Free Cryo-EM Gating Kinematics Benchmark (PRED-0012)."""
from __future__ import annotations

import ast
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cryoem_gating_kinematics_benchmark as bench  # noqa: E402
import independent_audit as audit  # noqa: E402


class TestCryoEMGatingKinematicsBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = bench.run_benchmark(run_falsify=True)
        cls.audit_report = audit.run_independent_audit()

    def test_ast_zero_shared_imports_in_independent_audit(self) -> None:
        src = (HERE / "independent_audit.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        forbidden = {"cryoem_gating_kinematics_benchmark", "tools", "subprocess", "re", "shlex"}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    root = alias.name.split(".")[0]
                    self.assertNotIn(root, forbidden)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    root = node.module.split(".")[0]
                    self.assertNotIn(root, forbidden)
                self.assertEqual(node.level, 0, "Relative imports forbidden in independent_audit.py")

    def test_structure_counts_and_construct_integrity(self) -> None:
        pc = self.report["positive_controls"]
        self.assertEqual(pc["n_5uak_ca"], 1139)
        self.assertEqual(pc["n_6msm_ca"], 1181)
        self.assertEqual(pc["n_common_ca"], 1130)
        self.assertEqual(self.audit_report["n_common_ca"], 1130)

    def test_nbd_dimer_closure_and_scaffolding_controls(self) -> None:
        pc = self.report["positive_controls"]
        self.assertGreater(pc["dimer_closure_pairs"]["464_1348"]["contraction_a"], 15.0)
        self.assertLess(pc["dimer_closure_pairs"]["464_1348"]["d_6msm"], 12.0)
        self.assertGreater(pc["dimer_closure_pairs"]["582_1246"]["contraction_a"], 15.0)
        self.assertLess(pc["dimer_closure_pairs"]["582_1246"]["d_6msm"], 12.0)
        for g_val in pc["gating_hinge_controls"].values():
            self.assertGreater(g_val, 0.45)
        for b_val in pc["rigid_beta_controls"].values():
            self.assertLess(b_val, 0.25)

    def test_se3_x_se3_invariance_and_swap_symmetry(self) -> None:
        pc = self.report["positive_controls"]
        self.assertLess(pc["max_se3_x_se3_drift"], 1e-12)
        self.assertEqual(pc["max_swap_symmetry_diff"], 0.0)

    def test_method1_and_method2_exact_agreement_across_all_1130_residues(self) -> None:
        seq = bench.load_uniprot_sequence()
        s5 = bench.parse_chain_a_ca(bench.DEFAULT_CIF_5UAK, "5UAK", bench.EXPECTED_5UAK_SHA256, seq)
        s6 = bench.parse_chain_a_ca(bench.DEFAULT_CIF_6MSM, "6MSM", bench.EXPECTED_6MSM_SHA256, seq)
        skin_float, deg_float = bench.compute_symmetric_lddt_kinematics(s5, s6)

        ca5_dec = audit.load_decimal_ca(audit.CIF_5UAK, "5UAK", seq)
        ca6_dec = audit.load_decimal_ca(audit.CIF_6MSM, "6MSM", seq)
        skin_frac, deg_frac = audit.compute_exact_fraction_s_kin(ca5_dec, ca6_dec)

        self.assertEqual(set(skin_float.keys()), set(skin_frac.keys()))
        for res_id in skin_float:
            self.assertEqual(deg_float[res_id], deg_frac[res_id])
            vf = skin_float[res_id]
            vq = skin_frac[res_id]
            if vf is None:
                self.assertIsNone(vq)
            else:
                self.assertIsNotNone(vq)
                self.assertAlmostEqual(vf, float(vq), places=12)

    def test_pred_0012_test_split_metrics_and_parity(self) -> None:
        self.assertEqual(self.report["outcome"], "supported")
        self.assertEqual(self.audit_report["outcome"], "supported")

        m1_cond = self.report["test"]["missense_s_kin"]["conditional_auroc"]
        m2_cond = self.audit_report["test_missense_s_kin_float"]
        self.assertAlmostEqual(m1_cond, m2_cond, places=6)
        self.assertEqual(self.audit_report["test_missense_s_kin_fraction"], "29/30")

        m1_abst = self.report["test"]["all_gated_g_kin"]["abstention_aware_auroc"]
        m2_abst = self.audit_report["test_gated_g_kin_float"]
        self.assertAlmostEqual(m1_abst, m2_abst, places=6)
        self.assertEqual(self.audit_report["test_gated_g_kin_fraction"], "76/91")

    def test_all_seven_planted_negative_controls_pass(self) -> None:
        fc = self.report["falsification_controls"]
        expected_keys = {
            "missing_e1371q_seq_dif_rejected",
            "unk_chain_collision_rejected",
            "empty_neighborhood_numeric_fallback_rejected",
            "asymmetric_state_swap_rejected",
            "non_missense_kinematic_score_rejected",
            "post_hoc_sign_flip_rejected",
            "codon_leakage_or_manual_source_rejected",
        }
        self.assertEqual(set(fc.keys()), expected_keys)
        for k, passed in fc.items():
            self.assertTrue(passed, f"Negative control {k} did not pass")

    def test_non_binary_labels_rejected_in_both_methods(self) -> None:
        with self.assertRaises(bench.CryoEMGatingKinematicsError):
            bench.evaluate_mann_whitney_auroc([1, 2, 0], [0.9, 0.5, 0.1])
        with self.assertRaises(ValueError):
            audit.wilcoxon_rank_sum_auroc([True, 0], [audit.Fraction(1, 2), audit.Fraction(1, 4)])

    def test_spatial_grid_exact_cell_boundary_15a(self) -> None:
        # (-1e-17, 0, 0) in cell -1 and (14.999999, 0, 0) in cell 0 / (15.0 + 2e-10, 0, 0) in cell +1
        g = bench.SpatialGrid15A({1: (-1e-17, 0.0, 0.0), 10: (14.999998, 0.0, 0.0)})
        self.assertEqual(g.neighbors_within(1, radius=15.0), {10})
        g_knife = bench.SpatialGrid15A({1: (-1e-17, 0.0, 0.0), 10: (15.0, 0.0, 0.0)})
        with self.assertRaises(bench.CryoEMGatingKinematicsError):
            g_knife.neighbors_within(1, radius=15.0)

    def test_catalog_allowlist_requires_verified_top_level_claim_label(self) -> None:
        import tempfile

        yaml_text = (
            "sources:\n"
            "  - id: missing-label-source\n"
            "    access: fetch\n"
            "  - id: notes-spoof-source\n"
            "    access: api\n"
            "    notes: >-\n"
            "      claim_label: verified\n"
            "  - id: genuine-verified-source\n"
            "    access: api\n"
            "    claim_label: verified\n"
        )
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "catalog.yaml"
            p.write_text(yaml_text, encoding="utf-8")
            admitted = bench.load_admitted_automated_catalog_ids(p)
            self.assertEqual(admitted, frozenset({"genuine-verified-source"}))

    def test_class_gated_score_rejects_non_missense_score_leakage(self) -> None:
        f508 = {"variant_id": "CFTR-F508del", "variant_family": "inframe_indel", "hgvs_p_or_c": "p.Phe508del"}
        with self.assertRaises(bench.CryoEMGatingKinematicsError):
            bench.class_gated_score(f508, 0.42)


if __name__ == "__main__":
    unittest.main()
