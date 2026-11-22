#!/usr/bin/env python3
"""Unit tests for Modulator-Free Apo 5UAK Packing Compliance Benchmark (PRED-0013)."""
from __future__ import annotations

import copy
import sys
import unittest
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cryoem_packing_compliance_benchmark as m1  # noqa: E402
import independent_audit as m2  # noqa: E402


class TestCryoEMPackingComplianceBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report_m1 = m1.run_benchmark(run_falsify=True)
        cls.report_m2 = m2.run_independent_audit()

    def test_01_cryptographic_and_construct_integrity(self) -> None:
        self.assertEqual(self.report_m1["positive_controls"]["placed_ca_5uak"], 1139)
        self.assertEqual(self.report_m2["placed_ca_5uak"], 1139)

    def test_02_collinear_four_point_exact_control(self) -> None:
        chain4 = {k: (3.8 * (k - 1), 0.0, 0.0) for k in (1, 2, 3, 4)}
        res = m1.compute_noncovalent_compliance_map(chain4)
        self.assertAlmostEqual(res[1] or 0.0, 36.0 / 49.0, places=12)
        self.assertAlmostEqual(res[2] or 0.0, 4.0 / 5.0, places=12)
        self.assertAlmostEqual(res[3] or 0.0, 4.0 / 5.0, places=12)
        self.assertAlmostEqual(res[4] or 0.0, 36.0 / 49.0, places=12)

    def test_03_exposed_vs_buried_and_train_positive_controls(self) -> None:
        pc = self.report_m1["positive_controls"]
        self.assertGreater(pc["min_exposed_control"], 0.050)
        self.assertLess(pc["max_buried_control"], 0.048)
        self.assertGreater(pc["min_exposed_control"], pc["max_buried_control"])
        self.assertGreater(pc["train_g178_s_flex"], pc["train_g85_s_flex"])

    def test_04_se3_invariance_and_so3_rational_zero_drift(self) -> None:
        self.assertLess(self.report_m1["positive_controls"]["max_se3_drift"], 1e-12)
        self.assertEqual(self.report_m2["so3_q_drift"], "0/1")
        self.assertEqual(self.report_m2["noncovalent_pairs_verified"], 646957)
        improper_reflection = (
            (2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0),
            (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
            (2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0),
        )
        with self.assertRaises(m1.CryoEMPackingComplianceError):
            m1.verify_proper_rotation_matrix(improper_reflection)

    def test_05_positive_scale_invariance_of_auroc(self) -> None:
        struct = m1.parse_5uak_chain_a()
        raw_coords = {k: pt for k, (_, pt) in struct.ca_by_seq.items()}
        annotations = m1.json.loads(m1.DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
        test_vars = [v for v in m1.validate_benchmark_annotations(annotations) if v["split"] == "test"]
        for r0 in (0.5, 1.0, 3.8, 12.5):
            cmap = m1.compute_noncovalent_compliance_map(raw_coords, r0_sq=r0 * r0)
            ev = m1.evaluate_split(test_vars, cmap, "trikafta_v23_responsive")
            self.assertAlmostEqual(ev["conditional_missense_auroc_s_flex"], 0.700000, places=9)
            self.assertAlmostEqual(ev["gated_abstention_auroc_g_flex"], 49.0 / 66.0, places=6)

    def test_06_wt_507_508_indels_receive_none_s_flex(self) -> None:
        by_name = {v["legacy_spl_name"]: v for v in self.report_m1["test_split"]["variants"]}
        self.assertIsNone(by_name["I507del"]["s_flex"])
        self.assertIsNone(by_name["F508del"]["s_flex"])
        self.assertEqual(by_name["I507del"]["g_flex"], 0.5)
        self.assertEqual(by_name["F508del"]["g_flex"], 0.5)

    def test_07_all_seven_falsification_controls_pass(self) -> None:
        fc = self.report_m1["falsification_controls"]
        self.assertEqual(len(fc), 7)
        for name, verdict in fc.items():
            self.assertEqual(verdict, "PASS", msg=name)

    def test_08_split_leakage_and_heldout_suffix_guards(self) -> None:
        annotations = m1.json.loads(m1.DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
        bad_suffix = copy.deepcopy(annotations)
        bad_suffix["variants"][0]["variant_family"] = "missense_heldout"
        with self.assertRaises(m1.CryoEMPackingComplianceError):
            m1.validate_benchmark_annotations(bad_suffix)

        bad_gap = copy.deepcopy(annotations)
        bad_gap["variants"][0]["codon_position"] = 452  # NBD1 domain + close to 455
        with self.assertRaises(m1.CryoEMPackingComplianceError):
            m1.validate_benchmark_annotations(bad_gap)

    def test_09_manual_sources_blocked_fail_closed(self) -> None:
        perms = m1.load_catalog_permissions()
        for blocked_id in ("mavedb-dump-zenodo", "cftr-dms-mckee-2023", "cftr-dms-tedman-2025"):
            with self.assertRaises(PermissionError):
                m1.assert_catalog_source_admitted(blocked_id, perms)

    def test_10_ast_independence_and_method_agreement(self) -> None:
        m2.verify_ast_independence(HERE / "independent_audit.py")
        te1 = self.report_m1["test_split"]
        te2 = self.report_m2["test"]
        tr1 = self.report_m1["train_split"]
        tr2 = self.report_m2["train"]
        self.assertEqual(te2["auc_cond_s_flex_fraction"], "7/10")
        self.assertEqual(te2["auc_abst_g_flex_fraction"], "49/66")
        self.assertEqual(te2["auc_abst_g_am_fraction"], "137/198")
        self.assertEqual(tr2["auc_cond_s_flex_fraction"], "1/3")
        self.assertEqual(tr2["auc_abst_g_flex_fraction"], "5/6")
        self.assertAlmostEqual(te1["conditional_missense_auroc_s_flex"], float(Fraction(7, 10)), places=6)
        self.assertAlmostEqual(te1["gated_abstention_auroc_g_flex"], float(Fraction(49, 66)), places=6)
        self.assertAlmostEqual(te1["conditional_missense_auroc_s_am"], float(Fraction(8, 15)), places=6)
        self.assertAlmostEqual(te1["gated_abstention_auroc_g_am"], float(Fraction(137, 198)), places=6)
        self.assertAlmostEqual(te1["sequence_consequence_auroc"], float(Fraction(15, 22)), places=6)
        self.assertAlmostEqual(tr1["conditional_missense_auroc_s_flex"], float(Fraction(1, 3)), places=6)
        self.assertAlmostEqual(tr1["gated_abstention_auroc_g_flex"], float(Fraction(5, 6)), places=6)
        self.assertEqual(self.report_m1["outcome"], "refuted")
        with self.assertRaises(m1.CryoEMPackingComplianceError):
            m1.directional_auroc_abstention([True, False], [0.9, 0.1])  # type: ignore[list-item]
        with self.assertRaises(m1.CryoEMPackingComplianceError):
            m1.directional_auroc_abstention([1, 0], [float("nan"), 0.1])


if __name__ == "__main__":
    unittest.main()
