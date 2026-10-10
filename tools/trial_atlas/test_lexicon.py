"""Tests for lexicon.py and lexicon.json: every rule's own examples, then each class with SYNTHETIC outcome wording (positives and
negatives), the vague-wording rule, the composite flag, safety subtypes, time-frame buckets and exact spans. No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lexicon  # noqa: E402
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402

LEX = lexicon.load()


def tag(measure, description="", time_frame=""):
    return LEX.tag_entry({"entry_id": "NCT00000000:P1", "nct_id": "NCT00000000", "measure": measure, "description": description,
                          "time_frame": time_frame})


def classes(measure, description=""):
    return [t["class"] for t in tag(measure, description)["tags"]]


# (class, wording that must be tagged with it, wording that must not). All SYNTHETIC.
CASES = [
    ("fev1", ["Absolute change in ppFEV1 through Week 24", "Forced expiratory volume in one second (L)"],
     ["Change in lung function", "Improvement in respiratory function", "Change in FVC"]),
    ("spirometry_other", ["Change in forced vital capacity", "FEF25-75 percent predicted", "Spirometry"],
     ["Change in ppFEV1 measured by spirometry"]),
    ("lci", ["Change in lung clearance index (LCI2.5)", "Multiple breath washout"], ["Mucociliary clearance rate"]),
    ("lung_function_unspecified", ["Lung function", "Change in pulmonary function tests"],
     ["Change in FEV1 (lung function)", "Liver function"]),
    ("sweat_chloride", ["Change in sweat chloride concentration"], ["Night sweats"]),
    ("npd", ["Nasal potential difference response"], ["Nasal polyp score"]),
    ("icm", ["Intestinal current measurement in rectal biopsies"], ["Intestinal obstruction episodes"]),
    ("cftr_biomarker_other", ["CFTR function in organoids"], ["CFTR function measured by sweat chloride", "CFTR modulator dose"]),
    ("exacerbations", ["Rate of protocol-defined pulmonary exacerbations"], ["Exacerbating factors"]),
    ("nutrition_growth", ["Change in BMI z-score", "Weight-for-age"], ["Bacterial growth in culture"]),
    ("gi_pancreatic", ["Coefficient of fat absorption", "Faecal elastase"], ["Pancreatitis"]),
    ("glucose_cfrd", ["Change in HbA1c", "Incidence of CF-related diabetes"], ["Glucagon-like peptide levels"]),
    ("liver", ["Change in liver stiffness"], ["Hepatic impairment cohort exposure"]),
    ("sinus_upper_airway", ["Sinus opacification", "SNOT-22 score"], ["Sinus rhythm on ECG"]),
    ("cfqr", ["CFQ-R respiratory domain score"], ["Quality of life questionnaire"]),
    ("pro_named_other", ["EQ-5D-5L index"], ["Quality of life"]),
    ("pro_unnamed", ["Quality of life"], ["CFQ-R respiratory domain", "EQ-5D-5L quality of life"]),
    ("imaging", ["Chest CT score", "Lung MRI ventilation defect"], ["Image quality questionnaire"]),
    ("microbiology_sputum", ["Sputum Pseudomonas aeruginosa density", "MRSA eradication"], ["Cell culture assay"]),
    ("inflammation_markers", ["Sputum neutrophil elastase", "Serum CRP"], ["Neutrophil function assay"]),
    ("healthcare_use", ["Number of hospitalizations", "Days of intravenous antibiotic treatment"], ["Hospital-based recruitment"]),
    ("survival_transplant", ["Time to lung transplantation", "All-cause mortality"],
     ["Number of participants with adverse events leading to death"]),
    ("safety_tolerability", ["Number of participants with adverse events", "Safety", "Tolerability"], ["Safely completed visits"]),
    ("pharmacokinetics", ["Maximum plasma concentration (Cmax)", "Half-life (t1/2)"], ["Lung clearance index", "Airway clearance"]),
    ("feasibility_adherence", ["Adherence to nebulised therapy", "Feasibility of recruitment"], ["Adherent mucus plugs"]),
]


class RuleExamplesTest(unittest.TestCase):
    def test_every_rule_has_examples_and_they_hold(self):
        for r in LEX.rules + LEX.composite:
            with self.subTest(rule=r["id"]):
                ex = r.get("examples", {})
                self.assertTrue(ex.get("match") and ex.get("no_match"), "every rule needs match and no_match examples")
                for t in ex["match"]:
                    self.assertIsNotNone(r["re"].search(t), f"{r['id']} should match {t!r}")
                for t in ex["no_match"]:
                    self.assertIsNone(r["re"].search(t), f"{r['id']} should not match {t!r}")

    def test_not_stated_rule_examples(self):
        for rid, rx in (("NS-01", LEX.placeholder), ("NS-V1", LEX.vague_ns)):
            r = next(x for x in LEX.data["not_stated_rules"] if x["id"] == rid)
            for t in r["examples"]["match"]:
                self.assertTrue(rx.search(t), f"{rid} {t!r}")
            for t in r["examples"]["no_match"]:
                self.assertFalse(rx.search(t), f"{rid} {t!r}")

    def test_every_class_except_other_and_not_stated_has_a_rule(self):
        ruled = {r["class"] for r in LEX.rules}
        self.assertEqual(set(LEX.class_order) - ruled, {"other", "not_stated"})

    def test_the_lexicon_is_versioned_and_hashed_canonically(self):
        self.assertTrue(LEX.version)
        self.assertEqual(len(LEX.sha256), 64)
        self.assertEqual(LEX.sha256, lexicon.Lexicon(LEX.data).sha256)

    def test_a_bad_lexicon_is_refused(self):
        bad = {**LEX.data, "rules": LEX.data["rules"] + [{**LEX.data["rules"][0], "class": "nope"}]}
        with self.assertRaises(ValueError):
            lexicon.Lexicon(bad)
        dup = {**LEX.data, "rules": LEX.data["rules"] + [LEX.data["rules"][0]]}
        with self.assertRaises(ValueError):
            lexicon.Lexicon(dup)


class ClassTest(unittest.TestCase):
    def test_each_class_positive_and_negative(self):
        for cls, positives, negatives in CASES:
            for text in positives:
                with self.subTest(cls=cls, text=text):
                    self.assertIn(cls, classes(text))
            for text in negatives:
                with self.subTest(cls=cls, negative=text):
                    self.assertNotIn(cls, classes(text))

    def test_vague_wording_is_never_fev1(self):
        for text in ("Lung function", "Change in pulmonary function", "Efficacy", "Respiratory function testing"):
            with self.subTest(text=text):
                self.assertNotIn("fev1", classes(text))
        t = tag("Lung function")
        self.assertTrue(t["flags"]["vague"])
        self.assertEqual(t["tags"][0]["class"], "lung_function_unspecified")

    def test_vague_measure_with_a_stated_measure_in_the_description(self):
        t = tag("Lung function", "Percent predicted FEV1 at Week 24.")
        self.assertEqual([x["class"] for x in t["tags"]], ["fev1"])
        self.assertEqual(t["tags"][0]["field"], "description")
        self.assertFalse(t["flags"]["vague"])

    def test_efficacy_alone_is_not_stated_and_vague(self):
        t = tag("Efficacy")
        self.assertEqual(t["status"], "not_stated")
        self.assertEqual(t["tags"][0]["rule_id"], "NS-V1")
        self.assertTrue(t["flags"]["vague"])

    def test_placeholder_is_not_stated(self):
        for text in ("TBD", "N/A", "", "   "):
            with self.subTest(text=text):
                t = tag(text)
                self.assertEqual(t["status"], "not_stated")
                self.assertEqual(t["tags"][0]["rule_id"], "NS-01")

    def test_unclassified_remainder(self):
        t = tag("Time to return to school")
        self.assertEqual(t["status"], "unclassified")
        self.assertEqual(t["tags"], [])

    def test_description_is_read_only_when_the_measure_gives_no_class(self):
        self.assertEqual(classes("Change in FEV1", "Safety will also be monitored."), ["fev1"])
        self.assertEqual(classes("Primary endpoint", "Change in sweat chloride from baseline."), ["sweat_chloride"])

    def test_composite_flag_and_components(self):
        t = tag("Composite of time to first pulmonary exacerbation or hospitalization")
        self.assertTrue(t["flags"]["composite"])
        self.assertFalse(t["flags"]["multi_class"])
        self.assertEqual([x["class"] for x in t["tags"]], ["exacerbations", "healthcare_use"])
        self.assertEqual(t["composite_span"]["rule_id"], "COMP-01")

    def test_multi_class_without_composite_wording(self):
        t = tag("Change in FEV1 and sweat chloride")
        self.assertTrue(t["flags"]["multi_class"])
        self.assertFalse(t["flags"]["composite"])

    def test_safety_subtypes(self):
        self.assertEqual(tag("Number of participants with serious adverse events")["safety_subtypes"], ["serious_adverse_events"])
        self.assertEqual(tag("Adverse events and serious adverse events")["safety_subtypes"], ["adverse_events", "serious_adverse_events"])
        self.assertEqual(tag("Safety")["safety_subtypes"], ["unspecified"])
        self.assertEqual(tag("Safety as measured by vital signs")["safety_subtypes"], ["vital_signs_ecg"])
        self.assertEqual(tag("Discontinuation due to intolerance")["safety_subtypes"], ["tolerability"])
        self.assertEqual(tag("Clinical laboratory abnormalities")["safety_subtypes"], ["laboratory"])
        self.assertEqual(tag("Dose-limiting toxicities")["safety_subtypes"], ["dose_limiting"])
        self.assertEqual(tag("Change in FEV1")["safety_subtypes"], [])

    def test_every_span_is_exact(self):
        texts = [p for _, pos, _ in CASES for p in pos] + ["Composite of exacerbation or hospitalization", "Efficacy", "TBD"]
        for text in texts:
            entry = {"entry_id": "NCT00000000:P1", "nct_id": "NCT00000000", "measure": text, "description": "", "time_frame": ""}
            for t in LEX.tag_entry(entry)["tags"]:
                with self.subTest(text=text, rule=t["rule_id"]):
                    self.assertEqual(entry[t["field"]][t["start"]:t["end"]], t["text"])
                    self.assertEqual(t["by"], "rule")
                    self.assertTrue(t["rule_id"])


class TimeFrameTest(unittest.TestCase):
    def test_buckets(self):
        cases = {
            "Pre-dose and up to 24 hours post-dose": "up to 1 day",
            "Day 28": "over 1 day to 4 weeks",
            "Up to 72 hours": "over 1 day to 4 weeks",
            "From baseline through Week 24": "over 4 weeks to 6 months",
            "6 months": "over 4 weeks to 6 months",
            "48 weeks": "over 6 to 12 months",
            "Up to 1 year": "over 6 to 12 months",
            "Twelve months": "over 6 to 12 months",
            "2 years": "over 12 months",
            "Baseline, Day 1, Day 15 and Month 18": "over 12 months",
            "TBD": "unparseable",
            "End of study": "unparseable",
            "": "not stated",
        }
        for text, bucket in cases.items():
            with self.subTest(text=text):
                self.assertEqual(LEX.timeframe(text)["timeframe_bucket"], bucket)

    def test_the_longest_duration_decides(self):
        self.assertEqual(LEX.timeframe("Week 4 and week 24")["timeframe_days"], 168.0)


class SnapshotTaggingTest(unittest.TestCase):
    def setUp(self):
        sf.block_network(self)

    def test_tag_snapshot_covers_in_scope_entries_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_cf_snapshot(Path(tmp) / "s")
            s = snap.load(Path(tmp) / "s")
            tags = lexicon.tag_snapshot(s, LEX)
        ids = [e["entry_id"] for e in tags["entries"]]
        self.assertEqual(len(ids), 20)
        self.assertNotIn("NCT00000011:P1", ids)                  # planted non-CF: out of scope, never tagged
        self.assertIn("NCT00000015:P0", ids)                      # no primary outcome: one not-stated entry
        self.assertEqual(tags["snapshot_sha256"], s.digest)
        self.assertEqual(tags["lexicon_sha256"], LEX.sha256)
        status = [e["status"] for e in tags["entries"]]
        self.assertEqual((status.count("rule"), status.count("not_stated"), status.count("unclassified")), (17, 2, 1))

    def test_non_cf_control_has_no_cf_specific_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_noncf_snapshot(Path(tmp) / "n")
            entries = lexicon.tag_all_entries(snap.load(Path(tmp) / "n"), LEX)
        self.assertEqual(len(entries), 8)
        self.assertEqual([e["entry_id"] for e in entries if any(t["class"] in LEX.cf_specific for t in e["tags"])], [])


if __name__ == "__main__":
    unittest.main()
