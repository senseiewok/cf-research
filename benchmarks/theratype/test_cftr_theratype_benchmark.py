#!/usr/bin/env python3
"""Unit tests for the Open CFTR Theratype Benchmark (`T-0925` / `PRED-0010`)."""
from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cftr_theratype_benchmark import (  # noqa: E402
    ALLOWED_FAMILIES,
    CftrVariantCase,
    DEFAULT_CATALOG_PATH,
    DEFAULT_FIXTURE_PATH,
    MCDONALD_2024_TABLE1_POSITIVE_CONTROLS,
    McDonaldVusControl,
    OPEN_CFTR_THERATYPE_BENCHMARK,
    SPLIT_VERSION,
    TheratypeBenchmarkError,
    audit_label_provenance,
    audit_split_leakage,
    evaluate_falsification_controls,
    evaluate_lofo_splits,
    evaluate_model_coverage_and_auroc,
    evaluate_pred_0010,
    load_benchmark_fixture,
    load_mavedb_cc0_csv,
    main,
    map_codon_to_p13569_domain,
    verify_positive_controls,
)
from independent_audit import audit_independently, main as independent_audit_main  # noqa: E402


class TestCftrTheratypeBenchmark(unittest.TestCase):
    def test_mcdonald_2024_table1_positive_controls(self) -> None:
        res = verify_positive_controls()
        self.assertEqual(res["status"], "PASS")
        self.assertEqual(res["controls_verified"], 7)

        bad_ctrl = list(MCDONALD_2024_TABLE1_POSITIVE_CONTROLS)
        bad_ctrl[0] = McDonaldVusControl("R31L", "p.Arg31Leu", 31, 0.91, "P", 0.17, "benign")
        with self.assertRaises(TheratypeBenchmarkError):
            verify_positive_controls(bad_ctrl)

    def test_split_leakage_domain_buffer_and_provenance(self) -> None:
        self.assertEqual(SPLIT_VERSION, "2026-11-21.v2")
        audit = audit_split_leakage(OPEN_CFTR_THERATYPE_BENCHMARK)
        self.assertTrue(audit["leakage_audit_passed"])
        self.assertEqual(audit["train_count"], 10)
        self.assertEqual(audit["test_count"], 20)
        self.assertEqual(audit["codon_position_overlap"], [])
        self.assertEqual(audit["positive_control_codon_overlap"], [])
        self.assertGreater(audit["min_cross_split_codon_gap"], 5)
        self.assertEqual(audit["structural_domain_overlap"], [])
        self.assertEqual(audit["invalid_families"], [])
        self.assertEqual(audit["provenance_violations"], [])
        self.assertEqual(audit["heldout_novel_families_in_test"], ["inframe_indel"])

        # Out-of-range / NaN missense score and frameshift token must fail provenance check
        nan_case = list(OPEN_CFTR_THERATYPE_BENCHMARK) + [
            CftrVariantCase("POISON-NAN", "p.Gly1349Asp", 1349, "missense", "test", 1, float("nan"), structural_domain="NBD2_C_terminal")
        ]
        self.assertFalse(audit_split_leakage(nan_case)["leakage_audit_passed"])

        fs_case = list(OPEN_CFTR_THERATYPE_BENCHMARK) + [
            CftrVariantCase("POISON-FS", "p.Lys684fs", 1349, "missense", "test", 0, 0.95, structural_domain="NBD2_C_terminal")
        ]
        self.assertFalse(audit_split_leakage(fs_case)["leakage_audit_passed"])

    def test_p13569_domain_bounds(self) -> None:
        self.assertEqual(map_codon_to_p13569_domain(1), "Lasso_N_terminal")
        self.assertEqual(map_codon_to_p13569_domain(508), "NBD1")
        self.assertEqual(map_codon_to_p13569_domain(1303), "NBD2_C_terminal")
        with self.assertRaises(TheratypeBenchmarkError):
            map_codon_to_p13569_domain(0)
        with self.assertRaises(TheratypeBenchmarkError):
            map_codon_to_p13569_domain(1481)

    def test_fda_spl_ground_truth_invariants(self) -> None:
        by_id = {c.variant_id: c for c in OPEN_CFTR_THERATYPE_BENCHMARK}
        self.assertEqual(by_id["CFTR-N1303K"].modulator_responsive, 1)
        self.assertEqual(by_id["CFTR-N1303K"].alyftrek_v9_responsive, 1)
        self.assertEqual(by_id["CFTR-N1303K"].kalydeco_v36_responsive, 0)
        for vid in ("CFTR-V520F", "CFTR-R560T", "CFTR-R1066C"):
            self.assertEqual(by_id[vid].modulator_responsive, 0)
            self.assertEqual(by_id[vid].alyftrek_v9_responsive, 1)
        for vid in ("CFTR-2789+5G>A", "CFTR-3849+10kbC>T"):
            self.assertEqual(by_id[vid].modulator_responsive, 1)
            self.assertEqual(by_id[vid].alyftrek_v9_responsive, 0)
        for vid in ("CFTR-F508del", "CFTR-I507del", "CFTR-G542X", "CFTR-W1282X", "CFTR-3849+10kbC>T"):
            self.assertIsNone(by_id[vid].alphamissense_score)
            self.assertIsNone(by_id[vid].revel_score)

    def test_pred_0010_refutation_and_baseline_comparison(self) -> None:
        pred = evaluate_pred_0010(OPEN_CFTR_THERATYPE_BENCHMARK)
        self.assertEqual(pred["prediction_id"], "PRED-0010")
        self.assertEqual(pred["outcome"], "refuted")

        am = pred["alphamissense_trikafta_v23"]
        base = pred["sequence_consequence_baseline_trikafta_v23"]
        self.assertEqual(am["test_total"], 20)
        self.assertEqual(am["test_scored"], 11)
        self.assertAlmostEqual(am["test_coverage"], 0.55, places=4)
        self.assertEqual(base["test_coverage"], 1.0)
        self.assertGreater(base["abstention_aware_auroc"], am["abstention_aware_auroc"])

        lofo = evaluate_lofo_splits(OPEN_CFTR_THERATYPE_BENCHMARK, "alphamissense")
        self.assertEqual(set(lofo.keys()), set(ALLOWED_FAMILIES))
        for fam, fold in lofo.items():
            self.assertEqual(fold["family_overlap_with_train"], [])
            if fam != "missense":
                self.assertEqual(fold["coverage"], 0.0)

    def test_manual_catalog_gate_strict_allowlist_and_yaml_hygiene(self) -> None:
        for manual_id in ("mavedb-dump-zenodo", "cftr-dms-mckee-2023", "cftr-dms-tedman-2025"):
            with self.assertRaises(PermissionError):
                load_mavedb_cc0_csv(DEFAULT_FIXTURE_PATH, DEFAULT_CATALOG_PATH, source_id=manual_id)

        with tempfile.TemporaryDirectory() as tmp:
            csv_file = Path(tmp) / "mavedb_mock.csv"
            csv_file.write_text("hgvs_pro,score\np.Gly551Asp,-1.20\n", encoding="utf-8")

            # Hostile YAML 1: inline comment after `access: manual # pending` + wrapped notes containing `access: fetch`
            cat_hostile_1 = Path(tmp) / "cat_hostile_1.yaml"
            cat_hostile_1.write_text(
                "sources:\n"
                "  - id: mavedb-dump-zenodo\n"
                "    access: manual # pending maintainer review\n"
                "    claim_label: verified\n"
                "    notes: >-\n"
                "      access: fetch\n"
                "      claim_label: verified\n",
                encoding="utf-8",
            )
            with self.assertRaises(PermissionError):
                load_mavedb_cc0_csv(csv_file, cat_hostile_1, source_id="mavedb-dump-zenodo")

            # Hostile YAML 2: `access: fetch` without `claim_label: verified` must fail closed
            cat_hostile_2 = Path(tmp) / "cat_hostile_2.yaml"
            cat_hostile_2.write_text(
                "sources:\n  - id: mavedb-dump-zenodo\n    access: fetch\n",
                encoding="utf-8",
            )
            with self.assertRaises(PermissionError):
                load_mavedb_cc0_csv(csv_file, cat_hostile_2, source_id="mavedb-dump-zenodo")

            # Valid admitted catalog entry (`access: fetch` + `claim_label: verified`) succeeds
            cat_admitted = Path(tmp) / "cat_admitted.yaml"
            cat_admitted.write_text(
                "sources:\n  - id: mavedb-dump-zenodo\n    access: fetch # admitted\n    claim_label: verified\n",
                encoding="utf-8",
            )
            rows = load_mavedb_cc0_csv(csv_file, cat_admitted, source_id="mavedb-dump-zenodo")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["hgvs_pro"], "p.Gly551Asp")

    def test_six_falsification_controls_and_cli(self) -> None:
        falsify = evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 6)
        self.assertEqual(falsify["total_controls"], 6)
        self.assertEqual(main(["--check", "--falsify"]), 0)
        self.assertEqual(independent_audit_main(["--check"]), 0)

    def test_independent_audit_exact_parity(self) -> None:
        ind = audit_independently(DEFAULT_FIXTURE_PATH)
        leak = audit_split_leakage(OPEN_CFTR_THERATYPE_BENCHMARK)
        am_eti = evaluate_model_coverage_and_auroc(OPEN_CFTR_THERATYPE_BENCHMARK, "alphamissense", "trikafta_v23")
        am_vtd = evaluate_model_coverage_and_auroc(OPEN_CFTR_THERATYPE_BENCHMARK, "alphamissense", "alyftrek_v9")
        base_eti = evaluate_model_coverage_and_auroc(
            OPEN_CFTR_THERATYPE_BENCHMARK, "sequence_consequence_baseline", "trikafta_v23"
        )
        base_vtd = evaluate_model_coverage_and_auroc(
            OPEN_CFTR_THERATYPE_BENCHMARK, "sequence_consequence_baseline", "alyftrek_v9"
        )

        self.assertEqual(ind["train_count"], leak["train_count"])
        self.assertEqual(ind["test_count"], leak["test_count"])
        self.assertEqual(ind["min_cross_split_codon_gap"], leak["min_cross_split_codon_gap"])
        self.assertEqual(ind["train_domains"], leak["train_domains"])
        self.assertEqual(ind["test_domains"], leak["test_domains"])

        for left, right in ((ind["alphamissense_trikafta_v23"], am_eti), (ind["alphamissense_alyftrek_v9"], am_vtd)):
            self.assertAlmostEqual(left["coverage_raw"], right["test_coverage_raw"], places=9)
            self.assertAlmostEqual(left["conditional_auroc_raw"], right["conditional_auroc_raw"], places=9)
            self.assertAlmostEqual(left["abstention_aware_auroc_raw"], right["abstention_aware_auroc_raw"], places=9)

        for left, right in (
            (ind["sequence_consequence_baseline_trikafta_v23"], base_eti),
            (ind["sequence_consequence_baseline_alyftrek_v9"], base_vtd),
        ):
            self.assertAlmostEqual(left["coverage_raw"], right["test_coverage_raw"], places=9)
            self.assertAlmostEqual(left["conditional_auroc_raw"], right["conditional_auroc_raw"], places=9)
            self.assertAlmostEqual(left["abstention_aware_auroc_raw"], right["abstention_aware_auroc_raw"], places=9)

    def test_label_provenance_separates_spl_listing_from_experimental_response(self) -> None:
        import copy

        _, _, payload = load_benchmark_fixture(DEFAULT_FIXTURE_PATH)
        clean = audit_label_provenance(payload)
        self.assertTrue(clean["audit_passed"], msg=str(clean))
        self.assertEqual(clean["panel_exposure_status"], "exposed_development_reference")
        self.assertEqual(clean["jurisdiction"], "US FDA")
        self.assertEqual(clean["canonical_transcript"], "ENST00000003084")
        self.assertEqual(clean["total_variants"], 30)
        self.assertEqual(clean["spl_listing_coverage"], 1.0)
        self.assertEqual(clean["experimental_measured_count"], 0)
        self.assertEqual(clean["experimental_unknown_count"], 30)
        self.assertEqual(clean["experimental_measured_coverage"], 0.0)
        self.assertEqual(clean["unsupported_assertion_count"], 0)
        self.assertEqual(
            clean["cross_spl_discordant_variants"],
            ["2789+5G->A", "3849+10kbC->T", "R1066C", "R560T", "V520F"],
        )

        # Planted defect 1: claiming measured experimental nonresponse without evidence source must fail
        poison_unsourced = copy.deepcopy(payload)
        poison_unsourced["variants"][0]["label_status"]["experimental_rescue_status"] = (
            "measured_experimental_nonresponse"
        )
        res_unsourced = audit_label_provenance(poison_unsourced)
        self.assertFalse(res_unsourced["audit_passed"])
        self.assertEqual(res_unsourced["unsupported_assertion_count"], 1)

        # Planted defect 2: binary 0 paired with listed_on_fda_spl must fail parity
        poison_parity = copy.deepcopy(payload)
        for v in poison_parity["variants"]:
            if v["trikafta_v23_responsive"] == 0:
                v["label_status"]["trikafta_v23_spl"] = "listed_on_fda_spl"
                break
        res_parity = audit_label_provenance(poison_parity)
        self.assertFalse(res_parity["audit_passed"])
        self.assertTrue(any("spl_parity_mismatch" in err for err in res_parity["violations"]))

        # Planted defect 3: marking exposed 2026-11-21.v2 panel as untouched_holdout must fail
        poison_exposure = copy.deepcopy(payload)
        poison_exposure["provenance"]["panel_exposure_status"] = "untouched_holdout"
        res_exposure = audit_label_provenance(poison_exposure)
        self.assertFalse(res_exposure["audit_passed"])
        self.assertTrue(
            any("exposed_panel_must_not_be_marked_untouched" in err for err in res_exposure["violations"])
        )


if __name__ == "__main__":
    unittest.main()
