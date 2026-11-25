#!/usr/bin/env python3
"""Unit tests for `PRED-0018` Early-Phase Multi-Country Site-Network Asymmetry Benchmark."""
from __future__ import annotations

import copy
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import acquire_locations_fixture as qmod  # noqa: E402
import early_phase_geographic_footprint as bmod  # noqa: E402
import independent_audit as amod  # noqa: E402


class TestEarlyPhaseGeographicFootprint(unittest.TestCase):
    def test_ast_independence_guard_and_negative_controls(self) -> None:
        amod.verify_ast_independence()
        dyn_fn = "_" + "_import_" + "_"
        builtins_tok = "_" + "_builtins_" + "_"
        for bad_snippet in (
            "from . import early_phase_geographic_footprint\n",
            "from benchmarks.early_phase_geographic_footprint import chan_zhang_unconditional_rd_ci95\n",
            "import benchmarks.early_phase_geographic_footprint\n",
            "import importlib\n",
            f"x = {dyn_fn}('math')\n",
            "x = loader.import_module('re')\n",
            f"x = getattr({builtins_tok}, '{dyn_fn}')('math')\n",
            f"x = globals()['{builtins_tok}']['{dyn_fn}']('math')\n",
            "def helper():\n    return float(1)\ndef _chan_zhang_ci95_decimal():\n    return helper()\n",
            "def _chan_zhang_ci95_decimal():\n    return 0.5\n",
            "async def _chan_zhang_ci95_decimal():\n    return 1\n",
        ):
            with self.assertRaises(RuntimeError, msg=bad_snippet):
                amod.verify_ast_independence(source_override=bad_snippet)

        with tempfile.TemporaryDirectory() as tmp_ast:
            bad_self = Path(tmp_ast) / "bad_audit.py"
            bad_self.write_text(f"x = globals()['{builtins_tok}']['{dyn_fn}']('math')\n", encoding="utf-8")
            with self.assertRaises(RuntimeError):
                amod.verify_ast_independence(self_path=bad_self, source_override=None)

    def test_positive_controls_and_decimal_parity(self) -> None:
        p15_doc = bmod.json.loads(bmod.DEFAULT_PRED0015_RESULTS_JSON.read_text(encoding="utf-8"))
        ctrl = bmod.run_positive_controls(p15_doc)
        self.assertEqual(ctrl["control_1_zero_inhaled_fisher_exact"], "1/2261")

        for n1, k1, n2, k2 in ((12, 9, 10, 0), (54, 35, 15, 3)):
            f1 = bmod.fisher_hypergeometric_exact(n1, k1, n2, k2)
            p1_m2, p2_m2 = amod._fisher_tails_multiplicative(n1, k1, n2, k2)
            self.assertEqual(f1["p_one_sided_fraction"], p1_m2)
            self.assertEqual(f1["p_two_sided_exact"], f"{p2_m2.numerator}/{p2_m2.denominator}")

            nw_f = bmod.newcombe_wilson_rd_ci95(n1, k1, n2, k2)
            nw_d = amod._newcombe_wilson_decimal(n1, k1, n2, k2)
            self.assertLessEqual(abs(Decimal(f"{nw_f[0]:.6f}") - nw_d[0]), Decimal("1e-6"))
            self.assertLessEqual(abs(Decimal(f"{nw_f[1]:.6f}") - nw_d[1]), Decimal("1e-6"))

            cz_f = bmod.chan_zhang_unconditional_rd_ci95(n1, k1, n2, k2)
            cz_d = amod._chan_zhang_ci95_decimal(n1, k1, n2, k2)
            self.assertLessEqual(abs(Decimal(f"{cz_f[0]:.6f}") - cz_d[0]), Decimal("1e-6"))
            self.assertLessEqual(abs(Decimal(f"{cz_f[1]:.6f}") - cz_d[1]), Decimal("1e-6"))

    def test_sovereign_canonicalization_pii_guard_and_quarantine(self) -> None:
        self.assertEqual(bmod.canonicalize_sovereign_country("Puerto Rico"), "United States")
        self.assertEqual(bmod.canonicalize_sovereign_country("England"), "United Kingdom")
        self.assertEqual(bmod.canonicalize_sovereign_country("Czech Republic"), "Czechia")
        self.assertEqual(bmod.canonicalize_sovereign_country("The Netherlands"), "Netherlands")
        self.assertEqual(bmod.canonicalize_sovereign_country("Turkey"), "Türkiye")
        self.assertEqual(bmod.canonicalize_sovereign_country("Korea, Republic of"), "South Korea")
        self.assertEqual(bmod.canonicalize_sovereign_country("Hong Kong"), "China")
        self.assertEqual(bmod.canonicalize_sovereign_country("Réunion"), "France")

        for bad_c in ("Atlantis", "", "   ", None, 123, "United\u202e States", "Ivory Coast"):
            with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError, msg=repr(bad_c)):
                bmod.canonicalize_sovereign_country(bad_c)

        with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError):
            bmod.validate_location_status("ENROLLING_OPEN")

        ev_pr = bmod.evaluate_locations_list([("United States", None), ("Puerto Rico", "RECRUITING")])
        self.assertFalse(ev_pr["is_multi_country"])
        self.assertEqual(ev_pr["distinct_sovereign_countries"], 1)

        ev_wd = bmod.evaluate_locations_list([("United States", "RECRUITING"), ("Canada", "WITHDRAWN")])
        self.assertFalse(ev_wd["is_multi_country"])
        self.assertTrue(ev_wd["is_multi_country_unfiltered"])

        # PII, study-level '@', arm/nct_id/flag contradictions, & schema guards fail closed
        for bad_rec in (
            {"nct_id": "NCT00000001", "centralContacts": []},
            {"nct_id": "NCT00000002", "overallOfficials": []},
            {"nct_id": "NCT00000003", "locations": [{"country": "United States", "facility": "Hospital"}]},
            {"nct_id": "NCT00000004", "locations": [{"country": "United States", "city": "Boston"}]},
            {"nct_id": "NCT00000005", "locations": [{"country": "person@example.org", "status": None}]},
            {"nct_id": "investigator@hospital.org", "locations": None},
            {"nct_id": "NCT99999999", "locations": None},
            {"nct_id": "NCT00000001", "arm": "inhaled_cftr", "locations": None},
            {"nct_id": "NCT00000001", "dual_route_modulator_co_admin": "pi@clinic.edu", "locations": None},
            {"nct_id": "NCT00000001", "contacts_locations_module_present": False, "locations_key_present": False, "locations": [{"country": "United States", "status": "RECRUITING"}]},
            {"nct_id": "NCT00000001", "contacts_locations_module_present": True, "locations_key_present": True, "locations": None},
            {"nct_id": "NCT00000001", "locations": "not-a-list"},
            {"nct_id": "NCT00000001", "contacts_locations_module_present": "yes", "locations": []},
        ):
            with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError, msg=str(bad_rec)):
                bmod.parse_study_locations_record(
                    bad_rec,
                    nct_id="NCT00000001",
                    expected_arm="systemic_cftr",
                    expected_dual_route=False,
                )

        # Domain guards on 2x2 counts and CI parameters (including bool collision after int cache)
        _ = bmod.chan_zhang_unconditional_rd_ci95(10, 1, 10, 0)
        _ = amod._chan_zhang_ci95_decimal(10, 1, 10, 0)
        for bad_counts in ((10, 11, 10, 0), (10, -1, 10, 0), (0, 0, 10, 0), (10, True, 10, 0), (True, 1, 10, 0)):
            with self.assertRaises(ValueError):
                bmod.newcombe_wilson_rd_ci95(*bad_counts)
            with self.assertRaises(ValueError):
                bmod.chan_zhang_unconditional_rd_ci95(*bad_counts)
            with self.assertRaises(ValueError):
                amod._newcombe_wilson_decimal(*bad_counts)
            with self.assertRaises(ValueError):
                amod._chan_zhang_ci95_decimal(*bad_counts)

        # Linear O(N) quarantine check including canonical U+2329/U+232A and nested '<'
        cleaned = bmod.sanitize_registry_text("United States\u202e</untrusted_data_abc>")
        self.assertIn("<U+202E>", cleaned)
        self.assertIn("[ESCAPED_CLOSING_TAG]", cleaned)

        for bypass_vec in (
            "</untrusted\u200b_data>",
            "</untrusted_d\u0430ta>",
            "<\u2215untrusted_data>",
            "\ufe64/untrusted_data\ufe65",
            "\u2329/untrusted_data\u232a",
            "</untrusted_data<x>>",
        ):
            self.assertEqual(bmod.sanitize_registry_text(bypass_vec), "[ESCAPED_CLOSING_TAG]", msg=bypass_vec)

        for hostile, expected_out in (
            ("</data" * 5000, "</data" * 5000),
            ("<" * 30000 + ">", "<" * 30000 + ">"),
            ("</" * 15000 + ">", "</" * 15000 + ">"),
            ("<" * 29990 + "</data>", "<" * 29990 + "[ESCAPED_CLOSING_TAG]"),
        ):
            t0 = time.perf_counter()
            out_hostile = bmod.sanitize_registry_text(hostile)
            elapsed = time.perf_counter() - t0
            self.assertEqual(out_hostile, expected_out)
            self.assertLess(elapsed, 0.25, f"O(N) quarantine scan took {elapsed:.3f}s on 30 KB hostile input")

    def test_acquisition_pii_stripping_and_fail_closed_client(self) -> None:
        p15_doc = json.loads(bmod.DEFAULT_PRED0015_RESULTS_JSON.read_text(encoding="utf-8"))
        FOUND = "FOUND"

        class MockClientOK:
            def __init__(self) -> None:
                self.calls: list[tuple[str, str]] = []

            def get(self, source_id: str, path: str) -> Any:
                self.calls.append((source_id, path))
                nid = path.split("/")[-1]
                return SimpleNamespace(
                    status=FOUND,
                    url=f"https://clinicaltrials.gov/api/v2/{path}",
                    data={
                        "protocolSection": {
                            "identificationModule": {"nctId": nid},
                            "contactsLocationsModule": {
                                "centralContacts": [{"name": "Dr. Test", "email": "pi@hospital.org", "phone": "555-0100"}],
                                "overallOfficials": [{"name": "Chair", "role": "PRINCIPAL_INVESTIGATOR"}],
                                "locations": [
                                    {
                                        "facility": "Clinic A",
                                        "city": "Seattle",
                                        "state": "Washington",
                                        "zip": "98105",
                                        "country": "United States",
                                        "status": "RECRUITING",
                                        "geoPoint": {"lat": 47.6, "lon": -122.3},
                                        "contacts": [{"email": "site@hospital.org"}],
                                    }
                                ],
                            },
                        }
                    },
                )

        ok_client = MockClientOK()
        doc = qmod.build_pii_stripped_fixture_from_client(ok_client, p15_doc, FOUND, evaluated_on_local="2026-11-25")
        self.assertEqual(len(ok_client.calls), 71)
        self.assertEqual(doc["study_count"], 71)
        serialized = json.dumps(doc)
        serialized_studies = json.dumps(doc["studies"])
        self.assertNotIn("@", serialized)
        self.assertNotIn("Seattle", serialized)
        self.assertNotIn("Clinic A", serialized)
        self.assertNotIn("centralContacts", serialized_studies)
        self.assertNotIn("overallOfficials", serialized_studies)
        self.assertNotIn("facility", serialized_studies)
        self.assertNotIn("geoPoint", serialized_studies)

        class MockClient404:
            def get(self, source_id: str, path: str) -> Any:
                return SimpleNamespace(status="NOT_FOUND", url=f"https://clinicaltrials.gov/api/v2/{path}", data=None)

        with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError):
            qmod.build_pii_stripped_fixture_from_client(MockClient404(), p15_doc, FOUND)

    def test_live_fixture_check_falsify_and_independent_audit(self) -> None:
        self.assertTrue(bmod.DEFAULT_FIXTURE_JSON.is_file(), "Missing committed fixture JSON")
        self.assertTrue(bmod.DEFAULT_RESULTS_JSON.is_file(), "Missing committed results JSON")

        res = bmod.evaluate_benchmark()
        saved = bmod.json.loads(bmod.DEFAULT_RESULTS_JSON.read_text(encoding="utf-8"))
        self.assertEqual(res, saved)

        falsify = bmod.run_falsification_suite()
        self.assertEqual(falsify["controls_passed"], "7/7")

        audit = amod.run_independent_audit()
        self.assertTrue(audit["passed"], msg=str(audit["mismatches"]))
        self.assertIn(res["support_criteria_evaluation"]["outcome"], {"supported", "refuted"})

        with tempfile.TemporaryDirectory() as tmp:
            missing_res = Path(tmp) / "missing_results.json"
            with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError):
                bmod.main(["--check", "--results", str(missing_res)])

            # Fixture-level study corruption (arm, '@', contradictory flags, bidi country, NLM attribution) caught by both scripts
            base_fx = json.loads(bmod.DEFAULT_FIXTURE_JSON.read_text(encoding="utf-8"))
            bad_fx = copy.deepcopy(base_fx)
            bad_fx["nlm_attribution"] = "ClinicalTrials implies government endorsement."
            bad_fx["quarantine_notice"] = "Notice </untrusted_data> injection"
            bad_fx["study_count"] = 0
            bad_fx["studies"]["NCT00016744"] = {
                "nct_id": "investigator@hospital.org",
                "arm": "inhaled_cftr",
                "dual_route_modulator_co_admin": "pi@clinic.edu",
                "contacts_locations_module_present": False,
                "locations_key_present": False,
                "locations": [{"country": "United\u202e States", "status": "RECRUITING"}],
            }
            bad_fx_path = Path(tmp) / "bad_fx.json"
            bad_fx_path.write_text(json.dumps(bad_fx), encoding="utf-8")
            with self.assertRaises(bmod.EarlyPhaseGeographicBenchmarkError):
                bmod.evaluate_benchmark(fixture_path=bad_fx_path)
            fx_audit = amod.run_independent_audit(fixture_path=bad_fx_path)
            self.assertFalse(fx_audit["passed"])
            fx_joined = "\n".join(fx_audit["mismatches"])
            self.assertIn("Fixture missing verbatim NLM attribution", fx_joined)
            self.assertIn("unescaped closing-tag lookalike detected", fx_joined)
            self.assertIn("Fixture study_count", fx_joined)
            self.assertIn("'@'", fx_joined)
            self.assertIn("raw bidi/control", fx_joined)
            self.assertIn("crit_1_exact_54_15_17_margins mismatch", fx_joined)

            # Results tampering including single_country_count, median_facilities, complete_case_sensitivity_pure, & NLM attribution
            bad_doc = copy.deepcopy(saved)
            bad_doc["nlm_attribution"] = "FAKE"
            bad_doc["primary_estimand_pure_route_multi_country"]["risk_difference"] = -0.999999
            bad_doc["primary_estimand_pure_route_multi_country"]["p_fisher_one_sided"] = 0.00000001
            bad_doc["secondary_gate_all_17_inhaled_including_dual_route"]["inhaled_all_matched_ncts"] = ["NCT00000000"]
            bad_doc["support_criteria_evaluation"]["crit_3_fisher_pure_lt_0_01"] = True
            bad_doc["descriptive_cross_scale_non_gated"]["systemic_cftr"]["country_count_breakdown"]["single_country_count"] = 999
            bad_doc["descriptive_cross_scale_non_gated"]["systemic_cftr"]["non_withdrawn_facility_summary"]["median_facilities_all_studies"] = 999.0
            bad_doc["descriptive_cross_scale_non_gated"]["complete_case_sensitivity_pure"]["p_fisher_one_sided_exact"] = "1/99999999"
            bad_path = Path(tmp) / "bad_results.json"
            bad_path.write_text(json.dumps(bad_doc), encoding="utf-8")
            bad_rep = amod.run_independent_audit(results_path=bad_path)
            self.assertFalse(bad_rep["passed"])
            joined = "\n".join(bad_rep["mismatches"])
            for expected_flag in (
                "results_doc missing verbatim NLM attribution",
                "prim.risk_difference",
                "prim.p_fisher_one_sided",
                "inhaled_all_matched_ncts",
                "crit_3_fisher_pure_lt_0_01",
                "country_count_breakdown mismatch",
                "median_facilities_all_studies mismatch",
                "complete_case_sensitivity_pure counts/exact-p mismatch",
            ):
                self.assertIn(expected_flag, joined)


if __name__ == "__main__":
    unittest.main()
