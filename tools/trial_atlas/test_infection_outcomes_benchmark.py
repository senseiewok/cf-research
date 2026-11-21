#!/usr/bin/env python3
"""Unit tests for tools/trial_atlas/infection_outcomes_benchmark.py & independent_audit.py (T-0944)."""
from __future__ import annotations

import importlib.util
import io
import json
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import infection_outcomes_benchmark as iob  # noqa: E402
import lexicon as lex  # noqa: E402
from care_team_trials import load_prisma_included_studies  # noqa: E402

IND_PATH = REPO_ROOT / "benchmarks" / "infection_outcomes" / "independent_audit.py"
_spec = importlib.util.spec_from_file_location("independent_audit", IND_PATH)
assert _spec and _spec.loader
ind_audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ind_audit)


class InfectionOutcomesBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = iob.run_infection_benchmark()
        raw_data = json.loads(iob.DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        cls.raw_by_nct = {str(s["nct_id"]): s for s in raw_data["studies"] if s.get("nct_id")}
        _, prisma = load_prisma_included_studies(iob.DEFAULT_SNAPSHOT_JSON)
        cls.prisma_ncts = {str(s["nct_id"]) for s in prisma if s.get("nct_id")}
        cls.lx = lex.Lexicon(json.loads(iob.DEFAULT_LEXICON_JSON.read_text(encoding="utf-8")))

    def test_snapshot_and_prisma_denominators_and_disjoint_positive_controls(self) -> None:
        self.assertEqual(self.report["snapshot_raw_interventional_count"], 1053)
        self.assertEqual(self.report["prisma_included_interventional_count"], 940)
        self.assertIn("does not evaluate safety or efficacy", self.report["nlm_disclaimer"])
        benchmark_ids = {gt["nct_id"] for gt in iob.INFECTION_BENCHMARK_20_GROUND_TRUTH}
        for pc_id in ("NCT02212587", "NCT01093521", "NCT00910351", "NCT01515007"):
            self.assertNotIn(pc_id, benchmark_ids)

    def test_strict_cohort_partitioning_and_taxonomy_accuracy(self) -> None:
        self.assertEqual(len(self.report["prisma_940_interventional_studies"]), 17)
        self.assertEqual(len(self.report["raw_1053_excluded_studies"]), 2)
        self.assertEqual(len(self.report["external_observational_studies"]), 1)

        by_prisma = {r["nct_id"]: r for r in self.report["prisma_940_interventional_studies"]}
        excl_by_id = {r["nct_id"]: r for r in self.report["raw_1053_excluded_studies"]}
        obs_by_id = {r["nct_id"]: r for r in self.report["external_observational_studies"]}

        # IGNITE (NCT02354859) and ELITE (NCT00391976) are P. aeruginosa antimicrobial trials, not NTM or diagnostic studies
        self.assertEqual(by_prisma["NCT02354859"]["pathogen_focus"], "pseudomonas_antimicrobial")
        self.assertEqual(by_prisma["NCT02354859"]["endpoint_locus"], "primary_clinical_fev1")
        self.assertEqual(by_prisma["NCT00391976"]["pathogen_focus"], "pseudomonas_antimicrobial")
        self.assertEqual(by_prisma["NCT00391976"]["endpoint_locus"], "primary_cfu_or_culture")

        # ABATE (NCT04294043) is X4-excluded NTM with secondary NTM clearance
        self.assertEqual(excl_by_id["NCT04294043"]["cohort_partition"], "raw_1053_x4_excluded")
        self.assertEqual(excl_by_id["NCT04294043"]["pathogen_focus"], "ntm_antimycobacterial")
        self.assertEqual(excl_by_id["NCT04294043"]["endpoint_locus"], "primary_safety_secondary_ntm_conversion")
        self.assertEqual(excl_by_id["NCT05616221"]["cohort_partition"], "raw_1053_x5_excluded")
        self.assertEqual(obs_by_id["NCT07312734"]["cohort_partition"], "external_observational")

    def test_out_of_sample_stop2_stewardship_and_modulator_microbiota_precedence(self) -> None:
        # Out-of-sample STOP2 (NCT02781610) must classify as antibiotic_stewardship without study-specific title overfitting
        stop2 = iob.classify_study_multi_field("NCT02781610", self.raw_by_nct, self.prisma_ncts, self.lx)
        self.assertEqual(stop2["cohort_partition"], "prisma_940_interventional")
        self.assertEqual(stop2["pathogen_focus"], "antibiotic_stewardship")
        self.assertEqual(stop2["endpoint_locus"], "primary_clinical_fev1")

        # Out-of-sample ETI intestine-lung microbiota study (NCT05937815) must preserve primary_cfu_or_culture locus
        eti_micro = iob.classify_study_multi_field("NCT05937815", self.raw_by_nct, self.prisma_ncts, self.lx)
        self.assertEqual(eti_micro["endpoint_locus"], "primary_cfu_or_culture")
        self.assertNotEqual(eti_micro["pathogen_focus"], "modulator_negative_control")

        # Out-of-sample general diagnostic swab/sputum concordance trial without study-specific title words
        synth_id = "NCT90000001"
        synth_map = {
            **self.raw_by_nct,
            synth_id: {
                "nct_id": synth_id,
                "study_type": "INTERVENTIONAL",
                "conditions": ["Cystic Fibrosis"],
                "brief_title": "Diagnostic Concordance of Oropharyngeal Swabs and Induced Sputum in Cystic Fibrosis",
                "official_title": "",
                "keywords": [],
                "interventions": [{"name": "Induced sputum and oropharyngeal swab", "type": "DIAGNOSTIC_TEST"}],
                "primary_outcomes": [{"measure": "Bacterial culture concordance", "description": "Sensitivity and specificity", "time_frame": "Day 1"}],
            },
        }
        diag_row = iob.classify_study_multi_field(synth_id, synth_map, self.prisma_ncts | {synth_id}, self.lx)
        self.assertEqual(diag_row["pathogen_focus"], "non_sputum_diagnostic_sampling")
        self.assertEqual(diag_row["endpoint_locus"], "primary_cfu_or_culture")

    def test_phage_primary_cfu_vs_phase1_safety_qpcr_distinction(self) -> None:
        by_id = {r["nct_id"]: r for r in self.report["prisma_940_interventional_studies"]}
        self.assertEqual(by_id["NCT06998043"]["endpoint_locus"], "primary_cfu_or_culture")
        self.assertTrue(by_id["NCT06998043"]["primary_cfu_ml_or_g_reported"])

        self.assertEqual(by_id["NCT07275905"]["endpoint_locus"], "primary_safety_secondary_qpcr")
        self.assertFalse(by_id["NCT07275905"]["primary_cfu_ml_or_g_reported"])
        self.assertIn("oral_swab", by_id["NCT07275905"]["specimen_modalities"])

    def test_multi_field_classifier_beats_primary_only_baseline(self) -> None:
        metrics = self.report["benchmark_metrics"]
        self.assertEqual(metrics["total_studies"], 20)
        self.assertEqual(metrics["multi_field_exact_matches"], 20)
        self.assertEqual(metrics["multi_field_exact_accuracy"], 1.0)
        self.assertEqual(metrics["held_out_error_fraction"], 0.0)
        self.assertEqual(metrics["primary_only_baseline_exact_matches"], 5)
        self.assertEqual(metrics["prisma_primary_microbiology_count"], 9)

    def test_dayman_2026_sample_vs_participant_counts_and_orientation(self) -> None:
        orient = self.report["registry_orientation"]
        self.assertEqual(orient["ntm_median_first_infection_age_years"], 26)
        self.assertNotIn("cff-research-priorities", self.report["verified_source_catalog_ids"])

        lit_by_doi = {r["doi"]: r for r in self.report["literature_anchors"]}
        dayman = lit_by_doi["10.1099/jmm.0.002150"]
        self.assertEqual(dayman["participants_enrolled"], 11)
        self.assertEqual(dayman["respiratory_pathogen_positive_mask_samples"], 9)
        self.assertEqual(dayman["mycobacterium_abscessus_positive_n"], 3)
        self.assertEqual(dayman["pseudomonas_aeruginosa_positive_mask_samples"], 0)
        self.assertEqual(dayman["post_clearance_16s_increased_yield_n"], 5)
        self.assertNotIn("respiratory_pathogen_positive_participants", dayman)

    def test_quarantine_blocks_greek_fullwidth_and_bidi_split_closing_tags(self) -> None:
        q_greek = iob.quarantine_untrusted_text("Payload </untrust\u03b5d_d\u03b1t\u03b1_9999> end", nonce="9999888877776666")
        self.assertGreaterEqual(q_greek["boundary_lookalikes_escaped"], 2)
        self.assertNotIn("</untrust\u03b5d_d\u03b1t\u03b1_9999>", q_greek["caged_text"])
        self.assertEqual(q_greek["caged_text"].count("</untrusted_data_9999888877776666>"), 1)

        q_bidi = iob.quarantine_untrusted_text("Payload </untrusted\u200b_data_9999888877776666> end", nonce="9999888877776666")
        self.assertGreaterEqual(q_bidi["boundary_lookalikes_escaped"], 2)
        self.assertGreaterEqual(q_bidi["trojan_chars_escaped"], 1)
        self.assertEqual(q_bidi["caged_text"].count("</untrusted_data_9999888877776666>"), 1)
        self.assertIn("&lt;/untrusted<U+200B>_data_9999888877776666&gt;", q_bidi["caged_text"])

        q_fw = iob.quarantine_untrusted_text("Payload \uFF1C/untrusted_data_9999\uFF1E end", nonce="9999888877776666")
        self.assertGreaterEqual(q_fw["boundary_lookalikes_escaped"], 2)

    def test_independent_audit_regex_free_and_measure_first_parity(self) -> None:
        self.assertNotIn("re", ind_audit.__dict__)
        counts = ind_audit.audit_independently()
        self.assertEqual(counts["raw_interventional_count"], 1053)
        self.assertEqual(counts["prisma_included_count"], 940)
        self.assertEqual(counts["directory_cross_verified_count"], 940)
        self.assertEqual(counts["prisma_primary_microbiology_in_20_benchmark"], 9)

        # Measure-first parity on adversarial inputs raised in blind review
        fev1_with_culture_desc = [{"measure": "Change in ppFEV1", "description": "In subjects with positive bacterial culture at baseline"}]
        self.assertFalse(ind_audit.has_primary_microbiology_measure_first(fev1_with_culture_desc))

        pa_eradication = [{"measure": "Eradication of P. aeruginosa", "description": ""}]
        self.assertTrue(ind_audit.has_primary_microbiology_measure_first(pa_eradication))

        fev1_pseudo = [{"measure": "FEV1 in patients with chronic Pseudomonas infection", "description": ""}]
        self.assertFalse(ind_audit.has_primary_microbiology_measure_first(fev1_pseudo))

        ae_achromo = [{"measure": "Adverse events in subjects with chronic Achromobacter infection", "description": ""}]
        self.assertFalse(ind_audit.has_primary_microbiology_measure_first(ae_achromo))

        neg_excluding = {"study_type": "INTERVENTIONAL", "conditions": ["Bronchiectasis excluding cystic fibrosis"]}
        self.assertEqual(ind_audit.classify_conditions_by_word_window(neg_excluding), "X5")

    def test_verify_benchmark_report_enforces_quarantine_cages_on_all_rows(self) -> None:
        tampered = dict(self.report)
        bad_rows = [dict(r) for r in self.report["prisma_940_interventional_studies"]]
        bad_rows[0]["quarantined_external_outcomes"] = {"nonce": "x", "caged_text": "uncaged"}
        tampered["prisma_940_interventional_studies"] = bad_rows
        with self.assertRaises(iob.InfectionBenchmarkError):
            iob.verify_benchmark_report(tampered)

    def test_falsification_and_cli_modes(self) -> None:
        f_rep = iob.run_falsification_suite()
        self.assertTrue(f_rep["falsify_passed"])
        self.assertEqual(len(f_rep["negative_controls"]), 5)

        buf = io.StringIO()
        with redirect_stdout(buf):
            self.assertEqual(iob.main(["--check"]), 0)
        self.assertIn("PASS infection_outcomes_benchmark", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
