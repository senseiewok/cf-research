#!/usr/bin/env python3
"""Unit tests for the Cryo-EM 3D potentiation-vs-folding pocket benchmark (PRED-0011)."""
from __future__ import annotations

import ast
import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cryoem_potentiation_benchmark as bench  # noqa: E402
import independent_audit as audit  # noqa: E402


class TestCryoEMPotentiationBenchmark(unittest.TestCase):
    def test_ast_zero_shared_imports_or_subprocess_in_independent_audit(self) -> None:
        src = (HERE / "independent_audit.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        forbidden = {"cryoem_potentiation_benchmark", "re", "shlex", "subprocess", "tools"}
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imported.add(alias.name.split(".")[0])
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        self.assertEqual(imported & forbidden, set())

    def test_positive_controls_and_se3_invariance(self) -> None:
        res = bench.run_positive_controls()
        self.assertTrue(res["passed"])
        self.assertEqual(res["placed_chain_a_residues"], 1162)
        self.assertEqual(res["potentiation_heavy_atoms"], 93)
        self.assertEqual(res["folding_heavy_atoms"], 103)
        self.assertEqual(res["pocket_controls_verified"], 8)
        self.assertLess(res["max_se3_drift"], 1e-12)
        self.assertGreater(res["train_g178_s3d"], res["train_g85_s3d"])

    def test_all_six_falsification_negative_controls_rejected(self) -> None:
        falsify = bench.run_falsification_suite()
        self.assertEqual(len(falsify), 6)
        for key, rejected in falsify.items():
            self.assertTrue(rejected, f"Negative control {key} did not fail closed")

    def test_catalog_allowlist_rejects_manual_and_unknown_sources(self) -> None:
        annot = json.loads(bench.DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
        variants = annot["variants"]
        for disallowed_src in ("cffpr-adr-2024", "cff-research-priorities", "closed-publisher-supp"):
            with self.assertRaises(PermissionError):
                bench.validate_variants_and_sources(variants, catalog_sources=["openfda-label", disallowed_src])

    def test_non_missense_rejects_unit_and_zero_shot_scores(self) -> None:
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.gated_score_from_missense_unit(
                {"variant_id": "CFTR-F508del", "variant_family": "inframe_indel", "hgvs_p_or_c": "p.Phe508del"},
                0.42,
            )
        annot = json.loads(bench.DEFAULT_ANNOTATIONS_JSON.read_text(encoding="utf-8"))
        for bad_key in ("s_3d", "alphamissense_score", "revel_score", "evo2_zero_shot_score"):
            mutated = [dict(v) for v in annot["variants"]]
            for row in mutated:
                if row["legacy_spl_name"] == "F508del":
                    row[bad_key] = 0.55
                    break
            with self.assertRaises(bench.CryoEMPotentiationError):
                bench.validate_variants_and_sources(mutated)

    def test_sequence_verifier_rejects_zero_index_and_wrong_declared_mutation(self) -> None:
        seq = bench.load_uniprot_sequence(verify_hash=False)
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.verify_ca_against_uniprot(
                {0: ("LEU", (0.0, 0.0, 0.0))},
                seq,
                declared_mutations={1371: "GLN"},
                declared_deletions=frozenset({508}),
            )
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.verify_ca_against_uniprot(
                {1371: ("TRP", (0.0, 0.0, 0.0))},
                seq,
                declared_mutations={1371: "GLN"},
                declared_deletions=frozenset({508}),
            )
        with self.assertRaises(AssertionError):
            audit.verify_ca_against_uniprot_strict(
                {0: ("LEU", (Decimal("0"), Decimal("0"), Decimal("0")))},
                seq,
                declared_mutations={1371: "GLN"},
                declared_deletions={508},
            )
        with self.assertRaises(AssertionError):
            audit.verify_ca_against_uniprot_strict(
                {1371: ("TRP", (Decimal("0"), Decimal("0"), Decimal("0")))},
                seq,
                declared_mutations={1371: "GLN"},
                declared_deletions={508},
            )

    def test_auroc_rejects_non_binary_labels_and_unequal_lengths_in_both_methods(self) -> None:
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.compute_auroc_metrics([1, 0, 2], [0.9, 0.1, 0.5])
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.compute_auroc_metrics([True, False], [0.9, 0.1])
        with self.assertRaises(ValueError):
            audit.wilcoxon_rank_sum_auroc([1, 0, 2], [Decimal("0.9"), Decimal("0.1"), Decimal("0.5")])
        with self.assertRaises(ValueError):
            audit.wilcoxon_rank_sum_auroc([1, 0], [Decimal("0.9")])

    def test_method1_and_method2_exact_agreement(self) -> None:
        m1 = bench.evaluate_benchmark()
        m2 = audit.run_independent_audit()

        ev = m1["pred_0011_evaluation"]
        self.assertEqual(ev["outcome"], "supported")
        self.assertTrue(ev["criterion_a_missense_passed"])
        self.assertTrue(ev["criterion_b_gated_passed"])

        self.assertEqual(m2["exact_fractions"]["missense_conditional_auroc_s3d"], "14/15")
        self.assertEqual(m2["exact_fractions"]["missense_conditional_auroc_sam"], "2/3")
        self.assertEqual(m2["exact_fractions"]["all_abstention_auroc_g3d"], "75/91")
        self.assertEqual(m2["exact_fractions"]["all_abstention_auroc_gam"], "67/91")
        self.assertEqual(m2["exact_fractions"]["all_abstention_auroc_seq_baseline"], "62/91")

        self.assertAlmostEqual(
            ev["test_missense_conditional_auroc_s3d"],
            m2["values"]["missense_conditional_auroc_s3d"],
            places=6,
        )
        self.assertAlmostEqual(
            ev["test_all_abstention_auroc_g3d"],
            m2["values"]["all_abstention_auroc_g3d"],
            places=6,
        )
        self.assertAlmostEqual(
            ev["test_all_abstention_auroc_seq_baseline"],
            m2["values"]["all_abstention_auroc_seq_baseline"],
            places=6,
        )

    def test_train_domain_limitation_explicitly_preserved(self) -> None:
        m1 = bench.evaluate_benchmark()
        train_s3d = m1["splits"]["train"]["missense_only"]["s_3d"]["conditional_auroc"]
        self.assertAlmostEqual(train_s3d, 0.25, places=6)

    def test_voxel_grid_matches_brute_force_on_all_placed_residues(self) -> None:
        seq = bench.load_uniprot_sequence(verify_hash=False)
        st = bench.parse_8eiq_structure(uniprot_seq=seq, verify_hash=False)
        pot_pts = [a.xyz for a in st.potentiation_atoms]
        fold_pts = [a.xyz for a in st.folding_atoms]
        pot_grid = bench.VoxelGrid3D(pot_pts, cell_size=8.0)
        fold_grid = bench.VoxelGrid3D(fold_pts, cell_size=8.0)
        self.assertEqual(len(st.ca_by_seq), 1162)
        for seq_id, (_, ca) in st.ca_by_seq.items():
            dp_grid = pot_grid.min_distance(ca)
            dp_brute = min(
                ((ca[0] - q[0]) ** 2 + (ca[1] - q[1]) ** 2 + (ca[2] - q[2]) ** 2) ** 0.5
                for q in pot_pts
            )
            df_grid = fold_grid.min_distance(ca)
            df_brute = min(
                ((ca[0] - q[0]) ** 2 + (ca[1] - q[1]) ** 2 + (ca[2] - q[2]) ** 2) ** 0.5
                for q in fold_pts
            )
            self.assertAlmostEqual(dp_grid, dp_brute, places=12, msg=f"pot mismatch at {seq_id}")
            self.assertAlmostEqual(df_grid, df_brute, places=12, msg=f"fold mismatch at {seq_id}")

    def test_auroc_rejects_single_class_and_out_of_range_gated_scores(self) -> None:
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.compute_auroc_metrics([1, 1, 1], [0.2, 0.5, 0.8])
        with self.assertRaises(bench.CryoEMPotentiationError):
            bench.gated_score_from_missense_unit({"variant_family": "missense", "hgvs_p_or_c": "p.Gly551Asp"}, 1.2)


if __name__ == "__main__":
    unittest.main()
