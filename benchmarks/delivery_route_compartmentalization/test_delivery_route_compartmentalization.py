#!/usr/bin/env python3
"""Unit and planted-tamper tests for `PRED-0015` (`delivery_route_compartmentalization`)."""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import delivery_route_compartmentalization as drc  # noqa: E402
import independent_audit as ia  # noqa: E402


class TestDeliveryRouteCompartmentalization(unittest.TestCase):
    def test_hypergeometric_calibration_and_unit_sum(self) -> None:
        res = drc.fisher_hypergeometric_exact(12, 9, 10, 0)
        self.assertEqual(res["p_one_sided_fraction"], Fraction(11, 24871))
        self.assertEqual(res["p_one_sided_exact"], "1/2261")
        self.assertEqual(res["support_sum_exact"], "1/1")

    def test_investigational_route_precedence_and_non_cftr_exclusions(self) -> None:
        vx522_dual = {
            "nct_id": "NCT05668741",
            "brief_title": "A Phase 1/2 Study of VX-522 in Participants With Cystic Fibrosis (CF)",
            "official_title": "",
            "interventions": [
                {"name": "VX-522 mRNA therapy", "type": "DRUG"},
                {"name": "IVA", "type": "DRUG"},
            ],
            "keywords": [],
        }
        cls_vx522 = drc.classify_study_regimen(vx522_dual)
        self.assertEqual(cls_vx522["arm"], "inhaled_cftr")
        self.assertTrue(cls_vx522["dual_route_modulator_co_admin"])
        self.assertEqual(ia._classify_regimen_independent(vx522_dual), ("inhaled_cftr", True))

        # Swapped brief_title / official_title VX-371 + Orkambi must still be excluded (Round-1 Finding 3b).
        vx371_swapped = {
            "nct_id": "NCT02709109",
            "brief_title": "Phase 2 Study of VX-371 With Orkambi",
            "official_title": "A Study to Evaluate the Safety and Efficacy of VX-371 in Subjects With Cystic Fibrosis",
            "interventions": [
                {"name": "VX-371 + HS", "type": "DRUG"},
                {"name": "Orkambi", "type": "DRUG"},
            ],
            "keywords": [],
        }
        self.assertEqual(drc.classify_study_regimen(vx371_swapped)["arm"], "non_cftr_or_other")
        self.assertEqual(ia._classify_regimen_independent(vx371_swapped), ("non_cftr_or_other", False))

        # Liposomal tobramycin in F508del CFTR patients must NOT leak into inhaled_cftr (Round-1 Finding 3a).
        liposomal_tobra = {
            "nct_id": "NCT99990003",
            "brief_title": "Phase 2 Study of Liposomal Tobramycin Inhalation in Adults With F508del CFTR",
            "official_title": "",
            "interventions": [{"name": "Liposomal tobramycin", "type": "DRUG"}],
            "keywords": [],
        }
        self.assertEqual(drc.classify_study_regimen(liposomal_tobra)["arm"], "non_cftr_or_other")
        self.assertEqual(ia._classify_regimen_independent(liposomal_tobra), ("non_cftr_or_other", False))

    def test_unflattened_schema_rejected_fail_closed(self) -> None:
        unflattened = {
            "nct_id": "NCT01225211",
            "study_type": "INTERVENTIONAL",
            "conditions": ["Cystic Fibrosis"],
            "phases": ["PHASE2"],
            "interventions": [{"name": "Ivacaftor", "type": "DRUG"}],
            "brief_title": "Study of Ivacaftor",
            "official_title": "Study of Ivacaftor",
            "primary_outcomes": [{"PrimaryOutcomeMeasure": "Sweat chloride", "timeFrame": "Day 28"}],
        }
        with self.assertRaises(drc.CompartmentalizationBenchmarkError):
            drc.assert_flattened_study_schema(unflattened)
        with self.assertRaises(drc.CompartmentalizationBenchmarkError):
            drc.classify_primary_outcomes({"nct_id": "NCT01225211", "primary_outcomes": ["Change in sweat chloride"]})

    def test_diagnostic_sweat_chloride_contamination_blocked(self) -> None:
        study = {
            "nct_id": "NCT99990009",
            "primary_outcomes": [
                {"measure": "Eligibility sweat chloride >= 60 mmol/L for CF diagnosis", "description": "", "time_frame": "Day 1"},
                {"measure": "Change from baseline in ppFEV1", "description": "", "time_frame": "Week 4"},
            ],
        }
        out = drc.classify_primary_outcomes(study)
        self.assertFalse(out["has_primary_sweat_chloride"])
        self.assertTrue(out["has_primary_airway_or_safety"])

    def test_snapshot_evaluation_and_falsification_suite(self) -> None:
        rep = drc.evaluate_snapshot()
        prim = rep["primary_estimand_early_phase"]
        self.assertEqual(rep["dataset"]["raw_interventional_count"], 1053)
        self.assertEqual(rep["dataset"]["prisma_included_count"], 940)
        self.assertEqual((prim["n_sys"], prim["k_sys"], prim["n_inh"], prim["k_inh"]), (54, 7, 17, 0))
        self.assertEqual(prim["p_fisher_one_sided_exact"], "160272/1203521")
        self.assertEqual(prim["outcome"], "refuted")
        self.assertEqual(prim["ledger_standing"], "refuted")

        f_rep = drc.run_falsification_suite()
        self.assertTrue(f_rep["falsify_passed"])
        self.assertEqual(f_rep["controls_passed"], 5)

    def test_independent_audit_exact_parity(self) -> None:
        audit_rep = ia.run_independent_audit()
        self.assertTrue(audit_rep["passed"], msg=str(audit_rep["mismatches"]))
        self.assertEqual(audit_rep["p_fisher_one_sided_exact"], "160272/1203521")
        self.assertEqual(audit_rep["p_fisher_two_sided_exact"], "20212/109411")
        self.assertEqual(audit_rep["chan_zhang_unconditional_ci95"], [-0.07838, 0.252184])

    def test_independent_audit_catches_planted_tamper_including_chan_zhang(self) -> None:
        orig = json.loads(drc.DEFAULT_RESULTS_JSON.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            tampered_counts = json.loads(json.dumps(orig))
            tampered_counts["primary_estimand_early_phase"]["k_sys"] = 25
            bad_path = Path(tmp) / "tampered_results.json"
            bad_path.write_text(json.dumps(tampered_counts), encoding="utf-8")
            rep = ia.run_independent_audit(results_path=bad_path)
            self.assertFalse(rep["passed"])
            self.assertTrue(any("2x2 counts mismatch" in m for m in rep["mismatches"]))

            # Round-1 Finding 1 regression: tampering with chan_zhang_unconditional_ci95 must fail!
            tampered_cz = json.loads(json.dumps(orig))
            tampered_cz["primary_estimand_early_phase"]["chan_zhang_unconditional_ci95"] = [0.5, 0.9]
            bad_cz_path = Path(tmp) / "tampered_cz.json"
            bad_cz_path.write_text(json.dumps(tampered_cz), encoding="utf-8")
            rep_cz = ia.run_independent_audit(results_path=bad_cz_path)
            self.assertFalse(rep_cz["passed"])
            self.assertTrue(any("Chan-Zhang CI95 mismatch" in m for m in rep_cz["mismatches"]))

            missing_keys_path = Path(tmp) / "missing_keys.json"
            missing_keys_path.write_text(json.dumps({"done_sha256": orig["done_sha256"]}), encoding="utf-8")
            rep_missing = ia.run_independent_audit(results_path=missing_keys_path)
            self.assertFalse(rep_missing["passed"])
            self.assertTrue(any("2x2 counts mismatch" in m for m in rep_missing["mismatches"]))


if __name__ == "__main__":
    unittest.main()
