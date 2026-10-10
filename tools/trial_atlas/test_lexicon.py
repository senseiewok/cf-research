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


class FixRoundLexiconTest(unittest.TestCase):
    """Review findings 2, 3, 8, 9 and 10 (fix round after e7ea5ed). Each case failed before its fix."""

    def test_multi_word_rules_tolerate_any_whitespace(self):
        multi = {}
        for r in LEX.rules + LEX.composite:
            for t in r["examples"]["match"]:
                if " " in t:
                    multi.setdefault(r.get("class", "composite"), []).append((r, t))
        classes_with_multiword_rules = {r.get("class", "composite") for r in LEX.rules + LEX.composite
                                        if " " in r["pattern"].replace("[- ]", "").replace("[\\s-]", "")}
        self.assertTrue(classes_with_multiword_rules <= set(multi), "every class with a multi-word rule needs a multi-word example")
        for cls, items in multi.items():
            for r, t in items:
                for variant in (t.replace(" ", "  ", 1), t.replace(" ", "\n", 1)):
                    with self.subTest(cls=cls, rule=r["id"], text=variant):
                        m = r["re"].search(variant)
                        self.assertIsNotNone(m, f"{r['id']} should match {variant!r}")
                        self.assertEqual(variant[m.start():m.end()], m.group(0))

    def test_widen_spaces_shapes(self):
        w = lexicon.widen_spaces
        self.assertEqual(w("a b"), r"a\s+b")
        self.assertEqual(w("(?<!non )x"), r"(?<!non\s)x")              # a lookbehind keeps a fixed width
        self.assertEqual(w("a ?b"), r"a\s?b")
        self.assertEqual(w("a[- ]b"), r"a[-\s]+b")
        self.assertEqual(w("a[- ]?b"), r"a(?:[-\s]+)?b")
        self.assertEqual(w(r"a\ b"), r"a\ b")                           # an escaped space stays as written

    def test_whitespace_variants_tag_with_exact_spans(self):
        t = tag("Change in sweat\nchloride")
        self.assertEqual(t["tags"][0]["class"], "sweat_chloride")
        self.assertEqual("Change in sweat\nchloride"[t["tags"][0]["start"]:t["tags"][0]["end"]], t["tags"][0]["text"])
        self.assertIn("fev1", classes("forced  expiratory volume in one second"))

    def test_non_serious_adverse_events_are_never_serious(self):
        for text in ("Number of non-serious adverse events", "Nonserious adverse events", "non serious adverse events",
                     "Incidence of non-SAEs"):
            with self.subTest(text=text):
                self.assertEqual(tag(text)["safety_subtypes"], ["adverse_events"])

    def test_time_frame_continuations(self):
        cases = {"Weeks 4-24": 168, "Week 4, 8, 12, 16": 112, "Months 1-6": 182.625, "Months 1, 3, 6": 182.625,
                 "Weeks 0, 4, 8, 12": 84, "Years 1 to 5": 1826.25, "Days 1-28": 28, "Days 1 through 29": 29,
                 "Weeks 2 and 6": 42, "30 minutes": 30 / 1440, "1 d": 1, "Minute 45": 45 / 1440}
        for text, days in cases.items():
            with self.subTest(text=text):
                self.assertAlmostEqual(LEX.timeframe(text)["timeframe_days"], round(days, 3), places=3)
        self.assertEqual(LEX.timeframe("Day 28")["timeframe_bucket"], "over 1 day to 4 weeks")
        self.assertEqual(LEX.timeframe("Day 29")["timeframe_bucket"], "over 4 weeks to 6 months")   # the written number decides
        self.assertEqual(LEX.timeframe("30 minutes")["timeframe_bucket"], "up to 1 day")

    def test_false_positives_removed(self):
        none_of = {
            "nutrition_growth": ["Dose by weight", "Weight-based dosing", "Molecular weight of the compound", "Body weight-based dosing",
                                 "Height of the nebulizer stand"],
            "microbiology_sputum": ["Sputum neutrophil elastase", "Sputum volume", "Daily sputum production",
                                    "Ivacaftor concentration in sputum"],
            "safety_tolerability": ["Exercise tolerance", "Heat tolerance", "Cold tolerance"],
            "feasibility_adherence": ["Peak oxygen uptake"],
            "survival_transplant": ["FEV1 in lung transplant recipients"],
            "liver": ["Liver function tests", "Liver enzymes"],
            "exacerbations": ["Rate of moderate or severe COPD exacerbations", "Asthma exacerbations", "ABPA exacerbations",
                              "Exacerbations of asthma"],
            "imaging": ["Bone mineral density by DXA", "Fat-free mass by DEXA"],
        }
        for cls, texts in none_of.items():
            for text in texts:
                with self.subTest(cls=cls, text=text):
                    self.assertNotIn(cls, classes(text))
        self.assertFalse(tag("Composite score of symptoms")["flags"]["composite"])
        self.assertFalse(tag("Brody composite score on chest CT")["flags"]["composite"])
        self.assertTrue(tag("Composite endpoint of death or transplant")["flags"]["composite"])
        self.assertEqual(tag("Liver function tests")["safety_subtypes"], ["laboratory"])
        self.assertIn("nutrition_growth", classes("Fat-free mass by DEXA"))
        self.assertIn("survival_transplant", classes("Time to lung transplantation"))
        self.assertIn("exacerbations", classes("Rate of pulmonary exacerbations"))
        self.assertIn("exacerbations", classes("Number of protocol-defined PEx"))

    def test_description_fallback_keeps_the_earliest_class_and_no_generic_safety(self):
        t = tag("Primary endpoint", "Change in sweat chloride; adverse events and FEV1 are also recorded.")
        self.assertEqual([x["class"] for x in t["tags"]], ["sweat_chloride"])
        self.assertTrue(t["tags"][0]["from_description"])
        self.assertEqual(tag("Primary endpoint", "Safety will be monitored throughout.")["status"], "unclassified")
        self.assertEqual(tag("Primary endpoint", "Adverse events will be recorded.")["status"], "unclassified")
        self.assertFalse(tag("Change in FEV1")["tags"][0].get("from_description"))

    def test_under_matches_added(self):
        has = {
            "pharmacokinetics": ["Plasma ivacaftor concentration", "Serum tobramycin concentrations", "Concentration of tobramycin in plasma",
                                 "Trough levels"],
            "healthcare_use": ["Days on IV antibiotics", "Antibiotic use", "Days of intravenous antibiotics"],
            "gi_pancreatic": ["Fecal fat excretion"],
            "microbiology_sputum": ["Positive culture for Pseudomonas", "Time to first positive culture"],
            "fev1": ["Change in pFEV1", "forced expiratory volume in 1s", "forced expiratory volume in 1 seconds",
                     "forced expiratory volume in the first second"],
            "spirometry_other": ["FEV (forced expiratory volume)"],
            "cftr_biomarker_other": ["CFTR gene expression"],
            "inflammation_markers": ["IL-10 concentration", "IL-1beta"],
            "exercise_capacity": ["Six-minute walk distance", "Peak oxygen uptake (VO2 peak)", "Exercise capacity", "Exercise tolerance"],
        }
        for cls, texts in has.items():
            for text in texts:
                with self.subTest(cls=cls, text=text):
                    self.assertIn(cls, classes(text))
        self.assertNotIn("pharmacokinetics", classes("Serum concentration of vitamin D"))
        self.assertNotIn("fev1", classes("FEV (forced expiratory volume)"))
        for text in ("Change in ALT", "AST elevations", "Total bilirubin", "Liver enzymes"):
            with self.subTest(text=text):
                self.assertEqual(tag(text)["safety_subtypes"], ["laboratory"])

    def test_the_negative_control_watches_exacerbations(self):
        self.assertIn("exacerbations", LEX.control_classes)
        self.assertTrue(LEX.cf_specific <= LEX.control_classes)

    def test_exercise_capacity_is_marked_as_a_taxonomy_addition(self):
        cls = next(c for d in LEX.data["domains"] for c in d["classes"] if c["id"] == "exercise_capacity")
        self.assertIn("maintainer", cls.get("added", ""))


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
