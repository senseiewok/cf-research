"""Tests for lexicon.py and lexicon.json: every rule's own examples, then each class with SYNTHETIC outcome wording (positives and
negatives), the vague-wording rule, the composite flag, safety subtypes, time-frame buckets and exact spans. No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import json
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
    ("sinus_upper_airway", ["Sinus opacification", "Chronic rhinosinusitis", "SNOT-22 score"], ["Sinus rhythm on ECG", "Snot on a tissue"]),
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
    # classes added in 0.3.0-draft
    ("lung_mechanics", ["Change in airway resistance by oscillometry", "Residual volume"], ["Change in FVC", "Antibiotic resistance"]),
    ("oxygenation", ["Nocturnal oxygen saturation by pulse oximetry", "Daytime PaCO2"], ["Peak oxygen uptake", "Tissue oxygen saturation"]),
    ("bone_density", ["Lumbar spine bone mineral density by DXA", "Change in BMD"], ["Body composition by DXA", "Bone age"]),
    ("mucociliary_clearance", ["Change in mucociliary clearance"], ["Lung clearance index", "Airway clearance technique"]),
    ("sputum_properties", ["24-hour sputum weight", "Sputum viscoelasticity"], ["Sputum culture", "Sputum neutrophil elastase"]),
    ("molecular_expression", ["Gene expression in blood leukocytes", "DNA methylation"], ["CFTR gene expression", "Genotype"]),
    ("oxidative_stress", ["Urinary isoprostane levels", "Oxidative stress markers"], ["Inhaled glutathione and FEV1"]),
    ("physical_activity", ["Daily step count", "Physical activity by accelerometry"], ["Exercise capacity", "Enzyme activity"]),
    ("muscle_strength", ["Handgrip strength", "Maximal inspiratory pressure"], ["Muscle mass", "MIP-1alpha levels"]),
]

ADDED_IN_0_3 = {"lung_mechanics", "oxygenation", "bone_density", "mucociliary_clearance", "sputum_properties", "molecular_expression",
                "oxidative_stress", "physical_activity", "muscle_strength"}


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
    """Whitespace, non-serious events, time-frame continuations, false positives and under-matches."""

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
                for variant in (t.replace(" ", "  "), t.replace(" ", "\n"), t.replace(" ", " \t\n ")):
                    with self.subTest(cls=cls, rule=r["id"], text=variant):
                        hit = LEX.find(r["re"], variant)
                        self.assertIsNotNone(hit, f"{r['id']} should match {variant!r}")
                        start, end, text = hit
                        self.assertEqual(variant[start:end], text)             # an exact piece of the ORIGINAL text

    def test_whitespace_is_normalised_with_an_offset_map(self):
        norm, spans = lexicon.normalise_whitespace("a  b\n\tc ")
        self.assertEqual(norm, "a b c ")
        self.assertEqual(spans[1], (1, 3))                                      # the collapsed run maps back to both spaces
        self.assertEqual(spans[2], (3, 4))

    def test_a_double_space_does_not_defeat_a_negative(self):
        # spacing must not defeat a negative
        self.assertEqual(tag("Number of participants with non  serious adverse events")["safety_subtypes"], ["adverse_events"])
        self.assertEqual(tag("Number of participants with non\nserious adverse events")["safety_subtypes"], ["adverse_events"])
        self.assertEqual(tag("Number of participants with serious  adverse events")["safety_subtypes"], ["serious_adverse_events"])
        self.assertEqual(classes("Exercise  tolerance"), ["exercise_capacity"])
        self.assertNotIn("nutrition_growth", classes("Dose adjusted by  body weight"))
        self.assertNotIn("nutrition_growth", classes("Dose adjusted by\tbody weight"))
        self.assertNotIn("safety_tolerability", classes("Oral glucose  tolerance test"))

    def test_long_digit_runs_in_a_time_frame_are_fast(self):
        import time
        start = time.perf_counter()
        LEX._tf_after.search("1" * 20000)
        self.assertLess(time.perf_counter() - start, 1.0)
        start = time.perf_counter()
        self.assertEqual(LEX.timeframe("1" * 20000)["timeframe_bucket"], "unparseable")    # over the length cap
        self.assertLess(time.perf_counter() - start, 1.0)

    def test_risky_patterns_are_refused_at_load(self):
        for bad in (r"(a[bc]*)*", r"(?:x+)+", r"(\w+ )+y", r"(?x) a b", r"(?ix)a", r"lung\ function"):
            data = json.loads(json.dumps(LEX.data))
            data["rules"][0]["pattern"] = bad
            with self.subTest(pattern=bad), self.assertRaises(ValueError):
                lexicon.Lexicon(data)
        data = json.loads(json.dumps(LEX.data))
        data["composite_rules"][0]["pattern"] = r"(a+)+"
        with self.assertRaises(ValueError):
            lexicon.Lexicon(data)

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


class Lexicon030Test(unittest.TestCase):
    """0.3.0-draft: the over-match fixes (the old false positive is gone, the true positive stays), the double tags resolved, the
    classes added, and the versioning notes. All wording is SYNTHETIC."""

    def test_version(self):
        # 0.3.1-draft and 0.3.2-draft build on 0.3.0-draft; their own pins are in Lexicon031Test and Lexicon032Test
        self.assertEqual(LEX.version, "0.3.2-draft")
        self.assertIn("versioning", LEX.data)

    def test_insulin_like_growth_factor_is_not_insulin(self):
        self.assertNotIn("glucose_cfrd", classes("Serum insulin-like growth factor 1"))
        self.assertNotIn("glucose_cfrd", classes("Insulin like growth factor binding protein"))
        self.assertIn("glucose_cfrd", classes("Insulin secretion after a meal"))

    def test_diabetes_only_as_an_outcome(self):
        for text in ("Change in FEV1 in adults with diabetes", "Diabetes insipidus", "Lung function in people living with diabetes"):
            with self.subTest(text=text):
                self.assertNotIn("glucose_cfrd", classes(text))
        for text in ("Incidence of CF-related diabetes", "Time to diagnosis of diabetes", "Diabetes onset", "Conversion to cystic fibrosis related diabetes"):
            with self.subTest(text=text):
                self.assertIn("glucose_cfrd", classes(text))

    def test_uptake_retention_compliance_outside_feasibility(self):
        for text in ("FDG uptake in the lungs", "Glucose uptake by muscle", "Uptake of FDG on PET", "CO2 retention", "Carbon dioxide retention",
                     "Lung compliance", "Static compliance", "Dynamic compliance"):
            with self.subTest(text=text):
                self.assertNotIn("feasibility_adherence", classes(text))
        for text in ("Uptake of the intervention", "Participant retention at 12 months", "Compliance with treatment", "Participant retention"):
            with self.subTest(text=text):
                self.assertIn("feasibility_adherence", classes(text))
        self.assertIn("lung_mechanics", classes("Static compliance"))
        self.assertIn("oxygenation", classes("Carbon dioxide retention"))

    def test_an_organism_naming_the_population_is_not_microbiology(self):
        self.assertEqual(classes("Change in FEV1 in patients with chronic Pseudomonas infection"), ["fev1"])
        self.assertEqual(classes("Quality of life in people with P. aeruginosa"), ["pro_unnamed"])
        for text in ("Sputum Pseudomonas aeruginosa density", "Time to Pseudomonas recurrence", "Eradication of S. aureus", "Pseudomonas acquisition",
                     "Time to first Pseudomonas infection"):
            with self.subTest(text=text):
                self.assertIn("microbiology_sputum", classes(text))

    def test_faecal_calprotectin_is_gastrointestinal_and_inflammation(self):
        # 0.3.0-draft tagged faecal calprotectin as inflammation only; 0.3.2-draft tags both (GI-04, Lexicon032Test, proposal)
        self.assertEqual(classes("Fecal calprotectin"), ["gi_pancreatic", "inflammation_markers"])
        self.assertEqual(classes("Faecal calprotectin concentration"), ["gi_pancreatic", "inflammation_markers"])
        self.assertEqual(classes("Sputum calprotectin"), ["inflammation_markers"])
        self.assertIn("gi_pancreatic", classes("Faecal elastase"))

    def test_snot22_is_named_and_sinus(self):
        # 0.3.0-draft tagged SNOT-22 as a named instrument only; 0.3.1-draft tags it as both (Lexicon031Test, proposal)
        for text in ("SNOT-22 score", "Sino-Nasal Outcome Test 22", "Change in SNOT 22"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["sinus_upper_airway", "pro_named_other"])
        self.assertEqual(classes("Sinus symptoms"), ["sinus_upper_airway"])

    def test_the_fev1_fvc_ratio_is_other_spirometry(self):
        for text in ("FEV1/FVC ratio", "Change in FEV1 / FVC", "Tiffeneau index FEV1/VC", "FEV1 to FVC ratio"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["spirometry_other"])
        self.assertEqual(classes("Change in FEV1 and FEV1/FVC"), ["fev1", "spirometry_other"])
        self.assertEqual(classes("FEV1/FVC and ppFEV1"), ["fev1", "spirometry_other"])
        self.assertEqual(classes("Absolute change in ppFEV1"), ["fev1"])

    def test_event_free_survival_is_not_survival(self):
        for text in ("Exacerbation-free survival", "Event-free survival", "Infection free survival"):
            with self.subTest(text=text):
                self.assertNotIn("survival_transplant", classes(text))
        self.assertIn("exacerbations", classes("Pulmonary exacerbation-free survival"))
        self.assertIn("survival_transplant", classes("Overall survival"))
        self.assertIn("survival_transplant", classes("Transplant-free survival"))

    def test_liver_function_is_safety(self):
        for text in ("Liver function", "Change in liver function"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["safety_tolerability"])
                self.assertEqual(tag(text)["safety_subtypes"], ["laboratory"])
        self.assertEqual(classes("Liver stiffness"), ["liver"])

    def test_sputum_amount_is_not_microbiology(self):
        for text in ("Wet sputum weight", "Volume of expectorated sputum", "Sputum dry weight", "Sputum rheology"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["sputum_properties"])
        self.assertEqual(classes("Sputum culture"), ["microbiology_sputum"])
        self.assertEqual(classes("Burkholderia spp. density in expectorated sputum"), ["microbiology_sputum"])

    def test_oximetry_named_explicitly_is_oxygenation(self):
        # second round: explicit oximetry inside a safety outcome gets both tags (Lexicon030Round2Test); a measure alone gets one
        self.assertEqual(classes("Safety measured by pulse oximetry"), ["oxygenation", "safety_tolerability"])
        self.assertEqual(classes("Change in oxygen saturation"), ["oxygenation"])
        self.assertNotIn("oxygenation", classes("Gas exchange variables"))
        self.assertNotIn("oxygenation", classes("Tissue oxygen saturation"))

    def test_named_instruments_and_unnamed_ratings(self):
        self.assertEqual(classes("Anxiety", "Measured by the Beck Anxiety Inventory."), ["pro_named_other"])
        self.assertEqual(classes("Dyspnoea score", "Modified Borg scale, 0 to 10."), ["pro_named_other"])
        self.assertEqual(classes("Pain score"), ["pro_unnamed"])
        self.assertEqual(classes("Change of maximum VAS"), ["pro_unnamed"])
        self.assertEqual(tag("Anxiety")["status"], "unclassified")      # a bare word names no instrument and no class

    def test_new_measures_in_existing_classes(self):
        has = {
            "exercise_capacity": ["One-minute sit-to-stand test repetitions", "1-minute STS", "Ventilatory threshold", "Incremental shuttle test"],
            "pharmacokinetics": ["Terminal elimination rate constant", "Fraction unbound", "Imipenem clearance", "Clearance of the study drug"],
            "nutrition_growth": ["Fat mass", "Daily energy intake", "Resting energy expenditure"],
            "inflammation_markers": ["Fractional exhaled nitric oxide (FeNO)"],
            "npd": ["Nasal transepithelial potential difference", "Nasal TEPD"],
            "cftr_biomarker_other": ["Beta-adrenergic sweat rate"],
            "imaging": ["Ventilation defect percentage", "Xenon MRI"],
            "feasibility_adherence": ["Recruitment rate", "Patient attrition", "Intervention fidelity"],
            "gi_pancreatic": ["Total intestinal transit time", "Gastric emptying"],
            "glucose_cfrd": ["CGM time in range", "Stimulated C-peptide"],
            "liver": ["Hepatic elasticity"],
            "lci": ["Nitrogen washout"],
            "spirometry_other": ["Change in FEV0.5"],
        }
        for cls, texts in has.items():
            for text in texts:
                with self.subTest(cls=cls, text=text):
                    self.assertIn(cls, classes(text))
        # 'chloride transport' alone names no test: the description decides
        self.assertEqual(tag("Chloride transport response")["status"], "unclassified")
        self.assertEqual(classes("Chloride transport response", "Measured by nasal TEPD."), ["npd"])
        self.assertNotIn("pharmacokinetics", classes("Clearance of Pseudomonas"))
        self.assertNotIn("pharmacokinetics", classes("Clearance of radiolabelled particles"))
        self.assertEqual(classes("CFTR gene expression"), ["cftr_biomarker_other"])

    def test_safety_additions(self):
        cases = {"Change in serum creatinine": ["laboratory"], "Acute kidney injury": ["laboratory"], "Audiometry": ["hearing"],
                 "Hearing loss": ["hearing"], "Body temperature": ["vital_signs_ecg"], "Respiratory rate": ["vital_signs_ecg"],
                 "QRS interval": ["vital_signs_ecg"], "Physical examination": ["physical_examination"], "Serum chemistry": ["laboratory"]}
        for text, subtypes in cases.items():
            with self.subTest(text=text):
                t = tag(text)
                self.assertEqual([x["class"] for x in t["tags"]], ["safety_tolerability"])
                self.assertEqual(t["safety_subtypes"], subtypes)
        self.assertNotIn("safety_tolerability", classes("AUC in participants with impaired renal function"))
        self.assertEqual(tag("Safety: physical examination")["safety_subtypes"], ["physical_examination"])

    def test_added_classes_are_marked_for_a_clinician_reader(self):
        added = {c["id"]: c for d in LEX.data["domains"] for c in d["classes"] if "0.3.0-draft" in c.get("added", "")}
        self.assertEqual(set(added), ADDED_IN_0_3)
        for cid, c in added.items():
            with self.subTest(cls=cid):
                self.assertIn("clinician", c["added"])
                self.assertLessEqual(len(c["gloss"]), 160)

    def test_unlike_cftr_gi_and_symptom_concepts_stay_apart(self):
        # the guardrail: pancreatic measures, symptoms and direct CFTR function are never one class
        self.assertNotEqual(LEX.domain_of["gi_pancreatic"], LEX.domain_of["cftr_biomarker_other"])
        self.assertEqual(classes("Faecal elastase"), ["gi_pancreatic"])
        self.assertEqual(classes("Sweat chloride"), ["sweat_chloride"])
        self.assertEqual(classes("Symptom severity"), ["pro_unnamed"])

    def test_changed_rules_carry_a_version_note(self):
        # the newest note comes first (0.3.1-draft), the 0.3.0-draft note is kept after it
        for r in LEX.data["rules"]:
            if "changed" in r:
                with self.subTest(rule=r["id"]):
                    self.assertTrue(r["changed"].startswith(("0.3.2-draft:", "0.3.1-draft:", "0.3.0-draft:")))
        changed = {r["id"] for r in LEX.data["rules"] if "0.3.0-draft:" in r.get("changed", "")}
        self.assertEqual(changed, {"FEV1-01", "FEV1-02", "SPIRO-02", "LUNGFN-V1", "EXAC-01", "GI-01", "GLU-01", "LIVER-01", "SINUS-01",
                                   "MICRO-01", "MICRO-03", "SURV-01", "SAFE-UNSPEC", "PK-01", "FEAS-01"})


class Lexicon030Round2Test(unittest.TestCase):
    """0.3.0-draft, second round: fixes from a blind adversarial review and a clinician-role review. All wording is SYNTHETIC."""

    def assert_has(self, cases):
        for cls, texts in cases.items():
            for text in texts:
                with self.subTest(cls=cls, text=text):
                    self.assertIn(cls, classes(text))

    def assert_not(self, cases):
        for cls, texts in cases.items():
            for text in texts:
                with self.subTest(cls=cls, text=text):
                    self.assertNotIn(cls, classes(text))

    def test_diabetes_outcomes_with_modifiers(self):
        self.assert_has({"glucose_cfrd": ["Incidence of type 2 diabetes", "Time to diabetes", "Number of participants who develop diabetes",
                                          "Diabetes-free survival", "Risk of developing diabetes", "Incidence of cystic fibrosis-related diabetes",
                                          "Onset of new-onset diabetes", "Prevalence of type 1 diabetes"]})
        self.assert_not({"glucose_cfrd": ["FEV1 in participants with a diagnosis of diabetes", "FEV1 in adults with type 2 diabetes",
                                          "Glucose-6-phosphate dehydrogenase activity", "Diabetes insipidus"]})

    def test_drug_clearance_is_positive_only(self):
        self.assert_has({"pharmacokinetics": ["Colistin clearance", "Tobramycin clearance", "Plasma clearance", "Apparent clearance",
                                              "Total clearance", "Renal clearance of the study drug", "Clearance of the study drug",
                                              "Clearance of tobramycin in plasma", "Clearance of the compound following a single oral dose",
                                              "Ceftazidime clearance"]})
        self.assert_not({"pharmacokinetics": ["Clearance of Haemophilus influenzae", "Clearance of nontuberculous mycobacteria",
                                              "Clearance of Stenotrophomonas maltophilia", "Clearance of the infection", "Clearance of symptoms",
                                              "Clearance of a radiolabelled molecule"]})

    def test_absorptive_clearance_is_left_to_other(self):
        for measure, desc in (("Absorptive clearance rate", "Clearance of a radiolabelled small molecule through absorption."),
                              ("DTPA absorption rate", "Total DTPA clearance minus mucociliary clearance rate.")):
            with self.subTest(measure=measure):
                t = tag(measure, desc)
                self.assertNotIn("pharmacokinetics", [x["class"] for x in t["tags"]])
                self.assertNotIn("mucociliary_clearance", [x["class"] for x in t["tags"]])
        self.assertEqual(classes("Mucociliary clearance rate"), ["mucociliary_clearance"])

    def test_glutathione_as_treatment_is_not_a_marker(self):
        self.assert_not({"oxidative_stress": ["Inhaled reduced glutathione", "Glutathione (GSH) inhalation", "Glutathione nebulization",
                                              "Intravenous glutathione", "Nebulised glutathione and FEV1"]})
        self.assert_has({"oxidative_stress": ["Sputum glutathione concentration", "Plasma glutathione levels", "GSH/GSSG ratio",
                                              "Glutathione redox status", "Plasma F2-isoprostanes"]})

    def test_sweat_secretion_needs_a_stimulation_or_cftr_context(self):
        self.assert_not({"cftr_biomarker_other": ["Sweat rate during exercise in the heat", "Transepidermal water loss by evaporimetry",
                                                  "Sweat secretion", "Skin water loss"]})
        self.assert_has({"cftr_biomarker_other": ["Beta adrenergic sweat secretion", "Cholinergic-stimulated sweat rate",
                                                  "Methacholine-stimulated sweat secretion", "CFTR-dependent sweat rate"],
                         "molecular_expression": ["CFTR-dependent gene expression in nasal cells"]})

    def test_lung_mechanics_needs_lung_context(self):
        self.assert_not({"lung_mechanics": ["Residual volume in the nebulizer cup", "Gastric residual volume during tube feeding",
                                            "Tidal volume set on the ventilator", "Tidal volume", "Inductive reactance of the device coil",
                                            "Grade X10 magnification images", "Change in specific image-based airway resistance (siRaw)"]})
        self.assert_has({"lung_mechanics": ["Specific airway resistance", "Respiratory system reactance at 5 Hz", "RV/TLC ratio",
                                            "Change in R5-R20", "Lung compliance", "Thoracic compliance", "Compliance of the lung",
                                            "Residual volume", "R5Hz by oscillometry"]})
        self.assertEqual(LEX.labels["lung_mechanics"], "Airway resistance, compliance and lung volumes")

    def test_feasibility_positive_forms(self):
        self.assert_not({"feasibility_adherence": ["Compliance of the lung", "Thoracic compliance", "CO2 retention", "CO₂ retention",
                                                   "Drug retention in the lung", "Uptake of 18F-FDG in the lungs", "Radiolabel uptake in the lung",
                                                   "Leukocyte recruitment", "Recruitment of neutrophils to the airway", "Alveolar recruitment manoeuvre",
                                                   "Oxygen uptake efficiency slope"]})
        self.assert_has({"feasibility_adherence": ["Compliance with the airway clearance regimen", "Medication compliance",
                                                   "Retention of participants at 12 months", "Participant retention", "Uptake of the intervention",
                                                   "Treatment uptake", "Participant recruitment", "Recruitment rate", "Study enrolment",
                                                   "Number of participants enrolled per month", "Dropout rate", "Adherence to therapy"]})

    def test_new_safety_rules_never_decide_from_the_description(self):
        for rid in ("SAFE-PE", "SAFE-VITAL2", "SAFE-LAB3", "SAFE-KIDNEY", "SAFE-KIDNEY2", "SAFE-HEARING"):
            with self.subTest(rule=rid):
                self.assertIn(rid, LEX.fallback_excluded)
        self.assertEqual(tag("Primary endpoint", "A physical examination is done at every visit.")["status"], "unclassified")

    def test_fev1_ratio_forms_are_not_fev1(self):
        for text in ("Change in FEV1-to-FVC ratio", "FEV1%FVC", "FEV1 over FVC", "FEV1/FEV6 ratio", "FEV1 / FVC (%)", "FEV1:FVC",
                     "Ratio of FEV1 to forced vital capacity",
                     "Forced expiratory volume in one second to forced vital capacity ratio",
                     "Forced expiratory volume in 1 second/forced vital capacity"):
            with self.subTest(text=text):
                self.assertNotIn("fev1", classes(text))
                self.assertIn("spirometry_other", classes(text))
        for text in ("FEV1 at week 24 and FEV1/FVC", "FEV1 to Week 24", "Forced expiratory volume in 1 second (FEV1)"):
            with self.subTest(text=text):
                self.assertIn("fev1", classes(text))
        self.assertNotIn("spirometry_other", classes("FEV 6 weeks after start"))
        self.assertIn("spirometry_other", classes("FEV0,5 z-score"))
        self.assertNotIn("fev1", classes("Lung function"))                     # vague wording is still never FEV1

    def test_smaller_false_positives(self):
        self.assert_not({"sputum_properties": ["Blood viscoelastic testing", "Blood rheology", "Sufficient sputum quantity for a pellet"],
                         "npd": ["TEPD of cultured cells", "TEPD in primary cell lines"],
                         "muscle_strength": ["Motor evoked potential (MEP) amplitude", "MEP latency", "MIP 1 alpha"],
                         "gi_pancreatic": ["Number of bowel movement-related symptoms", "Bowel movement symptom score"],
                         "safety_tolerability": ["FEV1 in participants with decreased renal function",
                                                 "Pharmacokinetics in participants with mild renal function impairment",
                                                 "Change in temperature of the reconstituted solution"]})
        self.assert_has({"sputum_properties": ["Grams of sputum expectorated", "Sputum percent solids", "Mucus solids percentage",
                                               "Expectorated sputum (g)", "Sputum wet weight in grams"],
                         "muscle_strength": ["Maximal expiratory pressure", "MIP and MEP", "Sniff nasal inspiratory pressure"],
                         "gi_pancreatic": ["Frequency of bowel movements", "Number of spontaneous bowel movements per week"],
                         "safety_tolerability": ["Creatinine clearance", "Decline in renal function", "Renal function tests"],
                         "npd": ["Nasal TEPD"]})

    def test_unicode_is_folded_and_spans_stay_exact(self):
        for text, cls in (("SpO₂ during sleep", "oxygenation"), ("PaCO₂", "oxygenation"), ("FEV₁/FVC", "spirometry_other"),
                          ("Change in ppFEV₁", "fev1"), ("Change in R5–R20", "lung_mechanics")):
            with self.subTest(text=text):
                t = tag(text)
                self.assertIn(cls, [x["class"] for x in t["tags"]])
                for x in t["tags"]:
                    self.assertEqual(text[x["start"]:x["end"]], x["text"])
        self.assertNotIn("glucose_cfrd", classes("Insulin–like growth factor 1"))          # en dash
        self.assertNotIn("glucose_cfrd", classes("Insulin‑like growth factor 1"))     # non-breaking hyphen
        self.assertNotIn("fev1", classes("FEV₁/FVC"))
        t = tag("SpO₂  during\nsleep")
        self.assertEqual("SpO₂  during\nsleep"[t["tags"][0]["start"]:t["tags"][0]["end"]], t["tags"][0]["text"])
        norm, spans = lexicon.normalise_whitespace("a₂–b")
        self.assertEqual(norm, "a2-b")
        self.assertEqual(spans, [(0, 1), (1, 2), (2, 3), (3, 4)])
        self.assert_has({"bone_density": ["Change in aBMD", "Tibial vBMD"]})

    def test_every_pattern_is_fast_on_adversarial_text(self):
        import time
        n = 20000
        texts = ["a" * n, "abc " * (n // 4), "clearance of " * (n // 13), "FEV1/" * (n // 5), "1" * n, "sputum " * (n // 7),
                 "incidence of a " * (n // 15), "mycin" * (n // 5), "PA-" * (n // 3), "R5-" * (n // 3), "a-" * (n // 2),
                 "oxygen uptake " * (n // 14), "a " * (n // 2), "-" * n,
                 "Abc-" * (n // 4), "Abc Def " * (n // 8)]   # hyphenated and spaced runs of capitalised words (PRO-06 was quadratic on them)
        rules = [(r["id"], r["re"]) for r in LEX.rules + LEX.composite] + [("NS-01", LEX.placeholder), ("NS-V1", LEX.vague_ns)]
        rules += [(f["id"], f["re"]) for f in LEX.fallback_restrictions]
        rules += [(r["id"], r["sup_re"]) for r in LEX.rules if r.get("sup_re")]
        for rid, rx in rules:
            for text in texts:
                norm, _ = lexicon.normalise_whitespace(text)
                start = time.perf_counter()
                rx.search(norm)
                elapsed = time.perf_counter() - start
                with self.subTest(rule=rid, text=text[:12]):
                    self.assertLess(elapsed, 0.1)

    def test_pk01_is_bounded_and_keeps_its_matches(self):
        for text in ("Plasma ivacaftor concentration", "Serum tobramycin concentrations", "Concentration of tobramycin in plasma",
                     "Ivacaftor concentration in sputum"):
            with self.subTest(text=text):
                self.assertIn("pharmacokinetics", classes(text))
        self.assertIn("{1,40}", next(r["pattern"] for r in LEX.data["rules"] if r["id"] == "PK-01"))

    def test_single_entry_alternatives_are_gone(self):
        pk3 = next(r["pattern"] for r in LEX.data["rules"] if r["id"] == "PK-03")
        for stem in ("derocol", "olozane", "zidime", "DTPA", "indium", "In-"):
            with self.subTest(stem=stem):
                self.assertNotIn(stem, pk3)
        pro4 = next(r["pattern"] for r in LEX.data["rules"] if r["id"] == "PRO-04")
        for name in ("GSAS", "PEI"):
            with self.subTest(name=name):
                self.assertNotIn(name, pro4)
        self.assertNotIn("Raw", next(r["pattern"] for r in LEX.data["rules"] if r["id"] == "LMECH-02").replace("sRaw", ""))
        self.assertNotIn("expectorated sputum \\(", next(r["pattern"] for r in LEX.data["rules"] if r["id"] == "SPUT-01"))

    def test_a_bare_exacerbation_measure_is_not_decided_by_other_outcomes_in_the_description(self):
        t = tag("Exacerbation rate", "Secondary measures: quality of life, Pseudomonas status and lung function.")
        self.assertEqual(t["status"], "unclassified")
        self.assertEqual(classes("Exacerbations frequency", "Number of respiratory exacerbations per year."), ["exacerbations"])
        self.assertEqual(classes("Number of participants with cystic fibrosis (CF) exacerbation"), ["exacerbations"])
        self.assertEqual(classes("Quality of life", "Exacerbation rate is also recorded."), ["pro_unnamed"])

    def test_timed_up_and_go_is_not_exercise_capacity(self):
        self.assertNotIn("exercise_capacity", classes("Timed up and go test"))
        self.assertIn("exercise_capacity", classes("30-second chair stand test"))
        self.assertIn("exercise_capacity", classes("One-minute sit-to-stand test"))

    def test_explicit_oximetry_in_a_safety_outcome_gets_both_tags(self):
        self.assertEqual(classes("Safety measured by pulse oximetry"), ["oxygenation", "safety_tolerability"])
        self.assertEqual(classes("Vital signs and pulse oximetry"), ["oxygenation", "safety_tolerability"])
        self.assertEqual(classes("Safety: oxygen saturation"), ["oxygenation", "safety_tolerability"])
        self.assertEqual(classes("Number of participants with desaturation adverse events"), ["safety_tolerability"])
        self.assertEqual(classes("Bronchospasm defined as oxygen desaturation after inhalation"), ["safety_tolerability"])
        self.assertEqual(classes("Nocturnal desaturation"), ["oxygenation"])

    def test_clinical_labels_and_glosses(self):
        self.assertEqual(LEX.labels["oxygenation"], "Oxygen and carbon dioxide levels")
        self.assertEqual(LEX.labels["muscle_strength"], "Muscle strength and endurance")
        self.assertEqual(LEX.labels["microbiology_sputum"], "Microbiology")
        dom = {d["id"]: d for d in LEX.data["domains"]}
        self.assertEqual(dom["microbiology_sputum"]["label"], "Microbiology")
        self.assertEqual(dom["pulmonary"]["label"], "Lung function and gas exchange")
        self.assertEqual(LEX.domain_of["oxidative_stress"], "molecular")
        self.assertEqual(LEX.domain_of["exercise_capacity"], "exercise_capacity")
        self.assertEqual(LEX.data["safety_subtype_labels"]["hearing"], "Hearing and balance")
        self.assertEqual(set(LEX.data["safety_subtype_labels"]), set(LEX.data["safety_subtypes"]))
        gloss = {c["id"]: c["gloss"] for d in LEX.data["domains"] for c in d["classes"]}
        self.assertIn("type-2", gloss["inflammation_markers"])
        self.assertIn("over time", gloss["mucociliary_clearance"])
        self.assertNotIn("how fast", gloss["mucociliary_clearance"].lower())
        self.assertIn("questionnaire", gloss["physical_activity"])
        self.assertIn("tissue", gloss["molecular_expression"])
        self.assertIn("z-score", gloss["bone_density"])
        self.assertNotIn("osteopenia", gloss["bone_density"])
        self.assertIn("percent solids", gloss["sputum_properties"])
        self.assertIn("exhaled", gloss["oxygenation"])
        for cid, g in gloss.items():
            with self.subTest(cls=cid):
                self.assertLessEqual(len(g), 160)

    def test_new_mechanisms_are_checked_at_load(self):
        for key, bad in (("suppressed_by_pattern", r"(a+)+"), ("measure_pattern", r"(?x) a b")):
            data = json.loads(json.dumps(LEX.data))
            if key == "suppressed_by_pattern":
                data["rules"][0]["suppressed_by_pattern"] = bad
            else:
                data["matching"]["fallback_restrictions"][0]["measure_pattern"] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                lexicon.Lexicon(data)
        data = json.loads(json.dumps(LEX.data))
        data["matching"]["fallback_restrictions"][0]["classes"] = ["nope"]
        with self.assertRaises(ValueError):
            lexicon.Lexicon(data)

    def test_suppressed_by_pattern_is_checked_in_the_same_field(self):
        self.assertEqual(classes("Mucociliary clearance rate", "Tc-sulfur colloid clearance; DTPA absorption was also measured."),
                         ["mucociliary_clearance"])


class Lexicon031Test(unittest.TestCase):
    """0.3.1-draft: rule-level fixes from a labelled development set (proposals/2026-10-11-lexicon-0.3.1-proposal.md). One test per
    changed or added rule, each with a positive and a negative. All wording is SYNTHETIC (made-up instrument and drug-free phrases)."""

    CHANGED = {"SPIRO-01", "SPIRO-02", "ICM-01", "CFTRBIO-02", "GI-01", "LIVER-01", "SINUS-01", "PRO-02", "MICRO-01", "MICRO-03",
               "SPUT-01", "INFL-01", "INFL-02", "OXS-01", "HCU-01", "SURV-01", "SURV-04", "PK-01", "PK-02", "PK-03", "FEAS-01", "FEAS-02"}
    ADDED = {"SINUS-02", "PRO-06"}

    def test_version_and_change_notes(self):
        self.assertTrue(LEX.version.startswith("0.3."))           # 0.3.2-draft notes come first; the 0.3.1 notes are kept after them
        noted = {r["id"] for r in LEX.data["rules"] if "0.3.1-draft:" in r.get("changed", "")}
        self.assertEqual(noted, self.CHANGED)
        for rid in self.ADDED:
            r = next(x for x in LEX.data["rules"] if x["id"] == rid)
            with self.subTest(rule=rid):
                self.assertIn("Added in 0.3.1-draft", r["note"])
        # no class added, removed or renamed
        self.assertEqual(len(LEX.class_order), 37)

    def test_pro06_named_instrument(self):
        for text in ("Example Coping Questionnaire (ECQ-II)", "Change in the Family Routine Inventory", "Partners in Example Scale score"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["pro_named_other"])
        for text in ("Pain on a Visual Analogue Scale", "Numeric Rating Scale for breathlessness", "Change in CF Questionnaire score"):
            with self.subTest(negative=text):
                self.assertNotIn("pro_named_other", classes(text))
        # the CFQ-R's own domains are called scales; a named activity questionnaire is physical activity
        self.assertEqual(classes("CFQ-R Respiratory Symptoms Scale"), ["cfqr"])
        self.assertEqual(classes("Physical activity by the Example Activity Estimation Scale"), ["physical_activity"])

    def test_pro02_gives_way_to_a_named_instrument_and_reads_self_report(self):
        self.assertEqual(classes("Coping measured by the Example Coping Questionnaire"), ["pro_named_other"])
        self.assertEqual(classes("Coping style", "Measured by the Example Flexibility Questionnaire."), ["pro_named_other"])
        self.assertEqual(classes("Self-reported sleep quality"), ["pro_unnamed"])
        self.assertEqual(classes("Quality of life questionnaire"), ["pro_unnamed"])
        self.assertNotIn("pro_unnamed", classes("Self-reporting bias in diaries"))

    def test_sinus01_sinus02_snot_is_named_and_sinus(self):
        self.assertEqual(classes("Sino-Nasal Outcome Test 22 score"), ["sinus_upper_airway", "pro_named_other"])
        self.assertEqual(classes("SNOT-20 total"), ["sinus_upper_airway", "pro_named_other"])
        self.assertNotIn("sinus_upper_airway", classes("Snotty nose reported by parents"))
        self.assertNotIn("sinus_upper_airway", classes("Sinus rhythm on ECG"))

    def test_spiro01_named_flows_drop_the_vague_tag(self):
        for text in ("FEF75 in percent predicted", "Pulmonary function: MEF50", "Slow vital capacity as a measure of lung function"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["spirometry_other"])
        self.assertEqual(classes("Lung function"), ["lung_function_unspecified"])
        self.assertNotIn("spirometry_other", classes("Vital status at one year"))

    def test_spiro02_forced_expiratory_manoeuvres(self):
        self.assertEqual(classes("Lung function", "Forced expiratory manoeuvres at each visit."), ["spirometry_other"])
        self.assertNotIn("spirometry_other", classes("Expiratory muscle training manoeuvres"))

    def test_icm01_needs_rectal_or_intestinal_context(self):
        self.assertEqual(classes("Short-circuit current in rectal biopsies"), ["icm"])
        self.assertEqual(classes("Intestinal current measurement"), ["icm"])
        self.assertNotIn("icm", classes("Short-circuit current across cultured airway cells"))

    def test_cftrbio02_takes_the_bare_chamber_measurement(self):
        self.assertEqual(classes("Short-circuit current across cultured airway cells"), ["cftr_biomarker_other"])
        self.assertEqual(classes("Channel function", "Assessed in an Ussing chamber."), ["cftr_biomarker_other"])
        self.assertNotIn("cftr_biomarker_other", classes("Short-circuit current in rectal biopsies"))

    def test_gi01_pancreatic_insufficiency_only_as_an_outcome(self):
        self.assertIn("gi_pancreatic", classes("Coefficient of nitrogen absorption"))
        self.assertIn("gi_pancreatic", classes("Incidence of exocrine pancreatic insufficiency"))
        self.assertEqual(classes("Safety of a new enzyme in participants with exocrine pancreatic insufficiency"), ["safety_tolerability"])

    def test_liver01_positive_forms(self):
        for text in ("Elastography of the liver", "Liver fat fraction", "Liver disease score"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["liver"] if "Elastography" not in text else ["liver", "imaging"])
        self.assertNotIn("liver", classes("Choline kinetics", "Bile acid output by the liver."))
        self.assertNotIn("liver", classes("Cod liver oil intake"))

    def test_micro03_sputum_as_a_sample_gives_way(self):
        self.assertEqual(classes("Sputum samples"), ["microbiology_sputum"])
        for text, cls in (("Sputum TNFα", "inflammation_markers"), ("Sputum IL-1ß", "inflammation_markers"),
                          ("Peak sputum concentration", "pharmacokinetics"), ("Sputum glutathione concentration", "oxidative_stress"),
                          ("Sputum expectorated per session", "sputum_properties")):
            with self.subTest(text=text):
                self.assertEqual(classes(text), [cls])

    def test_micro01_colonisation_and_colony_counts(self):
        self.assertEqual(classes("Nasal bacterial colonisation"), ["microbiology_sputum"])
        self.assertEqual(classes("Sputum colony count after treatment"), ["microbiology_sputum"])
        self.assertNotIn("microbiology_sputum", classes("Colony stimulating factor levels"))

    def test_sput01_sputum_expectorated(self):
        self.assertEqual(classes("Sputum expectorated during the session"), ["sputum_properties"])
        self.assertNotIn("sputum_properties", classes("Sputum expectorated for culture"))

    def test_infl01_cell_counts_and_proteases(self):
        for text in ("Blood neutrophils", "Quantification of blood neutrophils", "PMN counts in nasal lavage", "Proteinase 3 activity",
                     "Cathepsin G in sputum", "Leukotriene B4 levels"):
            with self.subTest(text=text):
                self.assertIn("inflammation_markers", classes(text))
        for text in ("Phagocytosis by blood neutrophils", "Gene expression in blood neutrophils", "Leukotriene receptor antagonist use"):
            with self.subTest(negative=text):
                self.assertNotIn("inflammation_markers", classes(text))

    def test_infl02_tnf_alpha(self):
        self.assertEqual(classes("Serum TNF-alpha"), ["inflammation_markers"])
        self.assertNotIn("inflammation_markers", classes("TNFR2 shedding"))

    def test_oxs01_nitrotyrosine_spelling_and_nrf2(self):
        self.assertEqual(classes("Plasma 3-nitrotyrosin"), ["oxidative_stress"])
        self.assertEqual(classes("Nrf2 activity in blood cells"), ["oxidative_stress"])
        self.assertNotIn("oxidative_stress", classes("Nrf2 activator dose escalation"))

    def test_hcu01_counts_not_anchors(self):
        for text in ("Number of antibiotic courses", "Time to first IV antibiotic course", "Number of clinic visits", "Length of hospitalization"):
            with self.subTest(text=text):
                self.assertIn("healthcare_use", classes(text))
        for text in ("FEV1 when the antibiotic course ends", "Change in lung function during hospitalization",
                     "Score at the first clinic visit", "Number of pulmonary exacerbations requiring hospitalization"):
            with self.subTest(negative=text):
                self.assertNotIn("healthcare_use", classes(text))
        # the composite with hospitalisation as its own component keeps both classes
        self.assertEqual(classes("Composite of time to first pulmonary exacerbation or hospitalization"), ["exacerbations", "healthcare_use"])

    def test_surv01_cell_survival_is_not_survival(self):
        self.assertEqual(classes("Survival"), ["survival_transplant"])
        for text in ("Neutrophil survival/apoptosis ratio", "Bacterial survival in sputum", "Cell survival in culture"):
            with self.subTest(text=text):
                self.assertNotIn("survival_transplant", classes(text))

    def test_surv04_bare_graft_is_not_a_transplant_outcome(self):
        self.assertEqual(classes("Lung transplant"), ["survival_transplant"])
        self.assertIn("survival_transplant", classes("Graft survival at one year"))
        self.assertNotIn("survival_transplant", classes("Graft colonisation by bacteria"))
        self.assertNotIn("survival_transplant", classes("Primary graft function score"))

    def test_pk01_pk02_peak_concentration_and_glucose_curves(self):
        for text in ("Peak sputum concentration", "Concentration of tobramycin in the sputum", "Area under the concentration-time curve"):
            with self.subTest(text=text):
                self.assertIn("pharmacokinetics", classes(text))
        for text in ("Plasma glucose area under the curve", "Area under the curve for plasma glucose", "Glucose AUC during the test",
                     "Serum concentration of inflammatory markers"):
            with self.subTest(negative=text):
                self.assertNotIn("pharmacokinetics", classes(text))

    def test_pk03_antibiotic_concentrations(self):
        self.assertIn("pharmacokinetics", classes("Sputum antibiotic levels after nebulisation"))
        self.assertNotIn("pharmacokinetics", classes("Antibiotic courses per year"))

    def test_feas01_misspelling(self):
        self.assertEqual(classes("Feasability evaluation"), ["feasibility_adherence"])
        self.assertNotIn("feasibility_adherence", classes("Feasting days"))

    def test_feas02_study_discontinuation_not_use_errors(self):
        self.assertEqual(classes("Children who discontinued the trial early"), ["feasibility_adherence"])
        self.assertNotIn("feasibility_adherence", classes("Critical use errors with the inhaler"))
        self.assertNotIn("feasibility_adherence", classes("Participants who discontinued the study drug"))

    # ---- DU-01: a generic patient-reported measure whose description names the instrument

    def test_du01_description_names_the_instrument(self):
        cases = (("Quality of life", "Measured with the CFQ-R at each visit.", ["cfqr"]),
                 ("Satisfaction scores", "Rated on the Example Therapy Rating Questionnaire.", ["pro_named_other"]),
                 ("Health-related quality of life", "Scored with the EQ-5D-5L.", ["pro_named_other"]),
                 ("Symptom burden", "The SNOT-22 total score.", ["sinus_upper_airway", "pro_named_other"]))
        for measure, desc, expected in cases:
            with self.subTest(measure=measure):
                t = tag(measure, desc)
                self.assertEqual([x["class"] for x in t["tags"]], expected)
                for x in t["tags"]:
                    self.assertTrue(x["from_description"])
                    self.assertEqual(desc[x["start"]:x["end"]], x["text"])          # the span is exact, in the description
                    self.assertEqual(x["field"], "description")

    def test_du01_negative_no_instrument_in_the_description(self):
        for desc in ("", "Self-rated on a scale from 0 to 10.", "A short questionnaire about daily life.", "Sinus symptoms are also asked about.",
                     "Pain on a Visual Analogue Scale.", "FEV1 is also recorded."):
            with self.subTest(desc=desc):
                t = tag("Quality of life", desc)
                self.assertEqual([x["class"] for x in t["tags"]], ["pro_unnamed"])
                self.assertFalse(t["tags"][0]["from_description"])

    def test_du01_never_touches_a_measure_that_named_another_class(self):
        # the measure names a non-PRO class: the description is not read, even if it names an instrument
        self.assertEqual(classes("Change in FEV1", "Quality of life by the CFQ-R is also measured."), ["fev1"])
        self.assertEqual(classes("FEV1 and quality of life", "Quality of life by the CFQ-R."), ["fev1", "pro_unnamed"])
        self.assertEqual(classes("Daily step count and satisfaction", "Satisfaction on the EQ-5D-5L."), ["pro_unnamed", "physical_activity"])
        # a measure that already names an instrument keeps its own classes
        self.assertEqual(classes("EQ-5D-5L index", "The CFQ-R is also given."), ["pro_named_other"])
        # a measure with no class still uses the ordinary earliest-match fallback
        self.assertEqual(classes("Primary endpoint", "Change in sweat chloride; the CFQ-R is also given."), ["sweat_chloride"])

    def test_du01_is_declared_and_checked_at_load(self):
        u = next(x for x in LEX.data["matching"]["description_upgrades"] if x["id"] == "DU-01")
        self.assertEqual(u["measure_classes"], ["pro_unnamed"])
        self.assertEqual(u["replaces"], ["pro_unnamed"])
        self.assertEqual({LEX.data["rules"][[r["id"] for r in LEX.data["rules"]].index(rid)]["class"] for rid in u["description_rules"]},
                         {"cfqr", "pro_named_other", "sinus_upper_airway"})
        for key, bad in (("description_rules", ["NOPE-01"]), ("description_rules", []), ("measure_classes", ["nope"]),
                         ("replaces", ["fev1"])):
            data = json.loads(json.dumps(LEX.data))
            data["matching"]["description_upgrades"][0][key] = bad
            with self.subTest(key=key, bad=bad), self.assertRaises(ValueError):
                lexicon.Lexicon(data)


class Lexicon032Test(unittest.TestCase):
    """0.3.2-draft: a second round from a larger development set (proposals/2026-10-11-lexicon-0.3.2-proposal.md). One test per changed
    or added rule or mechanism, each with a positive and a negative. All wording is SYNTHETIC."""

    CHANGED = {"PK-01", "PK-02", "GI-03", "PRO-01", "PRO-04", "PRO-06", "INFL-01", "SPUT-01", "FEAS-01", "FEAS-02", "MICRO-03"}
    ADDED = {"PK-04", "PK-05", "MCC-02", "GI-04", "GI-05"}

    def test_version_and_change_notes(self):
        self.assertEqual(LEX.version, "0.3.2-draft")
        noted = {r["id"] for r in LEX.data["rules"] if r.get("changed", "").startswith("0.3.2-draft:")}
        self.assertEqual(noted, self.CHANGED)
        for rid in self.ADDED:
            r = next(x for x in LEX.data["rules"] if x["id"] == rid)
            with self.subTest(rule=rid):
                self.assertIn("Added in 0.3.2-draft", r["note"])
        self.assertEqual(len(LEX.class_order), 37)

    def test_word_repetitions_in_changed_and_added_rules_are_bounded(self):
        import re
        for r in LEX.data["rules"]:
            if r["id"] in self.CHANGED | self.ADDED:
                with self.subTest(rule=r["id"]):
                    self.assertIsNone(re.search(r"\[[^\]]*\][+*]", r["pattern"]), "a character class repeated without a bound")

    def test_pk01_pk04_pk05_area_under_the_curve_is_weak(self):
        for text in ("Area under the plasma concentration-time curve", "Area under the curve from time zero to 24 hours", "AUC0-24",
                     "Plasma area under the curve of the study drug"):
            with self.subTest(text=text):
                self.assertEqual(classes(text), ["pharmacokinetics"])
        for text, cls in (("Area under the curve of MCC over 90 minutes", "mucociliary_clearance"), ("FEV1 area under the curve 0-4 hours", "fev1"),
                          ("Small bowel water content AUC", "gi_pancreatic"), ("Glucose AUC during the tolerance test", "glucose_cfrd")):
            with self.subTest(negative=text):
                self.assertEqual(classes(text), [cls])

    def test_micro03_gives_way_to_a_weak_drug_curve(self):
        self.assertEqual(classes("Sputum AUC of the inhaled drug"), ["pharmacokinetics"])
        self.assertEqual(classes("Sputum samples"), ["microbiology_sputum"])

    def test_pk02_keeps_its_other_acronyms(self):
        self.assertEqual(classes("Cmax and Tmax"), ["pharmacokinetics"])
        self.assertNotIn("pharmacokinetics", classes("cmax"))

    def test_mcc02_abbreviation(self):
        self.assertEqual(classes("Whole-lung MCC at 60 minutes"), ["mucociliary_clearance"])
        self.assertNotIn("mucociliary_clearance", classes("mcc of the device"))

    def test_gi03_gut_water_content(self):
        self.assertEqual(classes("Colonic water content by MRI"), ["gi_pancreatic", "imaging"])
        self.assertNotIn("gi_pancreatic", classes("Water content of the nebuliser chamber"))

    def test_gi04_gi05_intestinal_inflammation_and_symptom_instruments(self):
        self.assertEqual(classes("Change in gut inflammation"), ["gi_pancreatic", "inflammation_markers"])
        self.assertEqual(classes("IBS Severity Scoring System total"), ["gi_pancreatic", "pro_named_other"])
        self.assertEqual(classes("Number of children with IBS"), ["gi_pancreatic"])
        self.assertEqual(classes("Reflux symptoms in infants"), ["gi_pancreatic"])
        for text in ("Airway inflammation", "Pulmonary calprotectin", "Reflux of the solution into the cup", "ibs"):
            with self.subTest(negative=text):
                self.assertNotIn("gi_pancreatic", classes(text))

    def test_pro01_pro04_spellings(self):
        self.assertEqual(classes("Peds QL generic score"), ["pro_named_other"])
        self.assertEqual(classes("Rome IV diagnostic questionnaire"), ["gi_pancreatic", "pro_named_other"])
        self.assertNotIn("pro_named_other", classes("Rome visit schedule"))

    def test_pro06_one_word_name_with_its_abbreviation(self):
        self.assertEqual(classes("Self-Example Questionnaire for Children (SEQ-C)"), ["pro_named_other"])
        self.assertEqual(classes("Coping Inventory (CI-2) total"), ["pro_named_other"])
        self.assertNotIn("pro_named_other", classes("Body Mass Index (BMI)"))
        self.assertNotIn("pro_named_other", classes("Likert Scale (LS) of comfort"))

    def test_infl01_sputum_elastase_and_cell_count(self):
        self.assertEqual(classes("Change in sputum elastase activity"), ["inflammation_markers"])
        self.assertEqual(classes("Neutrophil cell count in lavage"), ["inflammation_markers"])
        self.assertNotIn("inflammation_markers", classes("Faecal elastase"))

    def test_sput01_wet_and_dry_weight(self):
        self.assertEqual(classes("Sputum wet and dry weight"), ["sputum_properties"])
        self.assertNotIn("sputum_properties", classes("Sputum culture and weight of the child"))

    def test_feas01_feas02_usability_and_reasons_for_leaving(self):
        for text in ("Usability of the inhaler app", "Ease of use of the device", "Participants who discontinued the study due to adverse events"):
            with self.subTest(negative=text):
                self.assertNotIn("feasibility_adherence", classes(text))
        self.assertEqual(classes("Example App Usability Scale"), ["pro_named_other"])
        self.assertEqual(classes("Children who discontinued the trial early"), ["feasibility_adherence"])
        self.assertEqual(classes("Acceptability of the intervention"), ["feasibility_adherence"])

    def test_du02_safety_measure_with_lung_tests_named_in_the_description(self):
        t = tag("Tolerability measured by pulmonary function tests", "Spirometry: FEV1 and FVC before and after the dose.")
        self.assertEqual([x["class"] for x in t["tags"]], ["fev1", "spirometry_other", "safety_tolerability"])
        self.assertFalse(t["flags"]["vague"])
        self.assertTrue(all(x["from_description"] for x in t["tags"] if x["class"] != "safety_tolerability"))
        # no named lung measure in the description: the vague tag stays
        self.assertEqual(classes("Tolerability measured by pulmonary function tests", "Done at each visit."),
                         ["lung_function_unspecified", "safety_tolerability"])
        # a safety-only measure never gets a lung class this way
        self.assertEqual(classes("Number of participants with adverse events", "FEV1 is also measured."), ["safety_tolerability"])
        # a vague measure alone still uses the ordinary fallback
        self.assertEqual(classes("Lung function", "FEV1 percent predicted."), ["fev1"])

    def test_fc01_two_marker_families_in_one_description(self):
        t = tag("Breath sample analysis", "Levels of nitrotyrosine and neutrophil elastase in the sample.")
        self.assertEqual([x["class"] for x in t["tags"]], ["inflammation_markers", "oxidative_stress"])
        self.assertTrue(all(x["from_description"] for x in t["tags"]))
        # a class outside the group is not added
        self.assertEqual(classes("Breath sample analysis", "Nitrotyrosine levels; FEV1 is also recorded."), ["oxidative_stress"])

    def test_fc02_gi_measure_and_its_questionnaire_in_one_description(self):
        self.assertEqual(classes("Mean score", "Gastroesophageal reflux disease symptoms scored with the Example Reflux Symptom Questionnaire."),
                         ["gi_pancreatic", "pro_named_other"])
        self.assertEqual(classes("Mean score", "Reflux symptoms; the CFQ-R is also given."), ["gi_pancreatic"])

    def test_new_mechanisms_are_checked_at_load(self):
        for bad in (["nope", "fev1"], ["fev1"]):
            data = json.loads(json.dumps(LEX.data))
            data["matching"]["fallback_companions"][0]["classes"] = bad
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                lexicon.Lexicon(data)


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

    # ---- tag output, the classes listing and the class labels

    def test_tag_never_overwrites_its_output(self):
        import contextlib
        import io
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_cf_snapshot(Path(tmp) / "s")
            out = Path(tmp) / "tags.json"
            out.write_text("KEEP", encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(lexicon.main(["tag", str(Path(tmp) / "s"), "--out", str(out)]), 2)
            self.assertEqual(out.read_text(encoding="utf-8"), "KEEP")

    def test_the_not_stated_class_is_labelled_no_measure_named(self):
        # this class label and its unchanged domain label; every label waits for a clinician reader
        self.assertEqual(LEX.labels["not_stated"], "No measure named")
        domain = next(d for d in LEX.data["domains"] if d["id"] == "not_stated")
        self.assertEqual(domain["label"], "Not stated")

    def test_classes_lists_every_class_with_its_domain_and_gloss(self):
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(lexicon.main(["classes"]), 0)
        lines = out.getvalue().strip().splitlines()
        for d in LEX.data["domains"]:
            for c in d["classes"]:
                with self.subTest(cls=c["id"]):
                    self.assertTrue(c.get("gloss", "").strip(), "every class needs a gloss")
                    self.assertLessEqual(len(c["gloss"]), 160)
                    self.assertTrue(any(line.split()[:2] == [c["id"], d["id"]] and c["gloss"] in line for line in lines))

    def test_non_cf_control_has_no_cf_specific_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            sf.build_noncf_snapshot(Path(tmp) / "n")
            entries = lexicon.tag_all_entries(snap.load(Path(tmp) / "n"), LEX)
        self.assertEqual(len(entries), 8)
        self.assertEqual([e["entry_id"] for e in entries if any(t["class"] in LEX.cf_specific for t in e["tags"])], [])


class Canonical8TaxonomyAndSpanTest(unittest.TestCase):
    def setUp(self):
        sf.block_network(self)

    def test_all_lexicon_classes_map_to_canonical_8(self):
        for c in LEX.class_order:
            canonical = lexicon.to_canonical_8(c)
            self.assertIn(canonical, lexicon.CANONICAL_8_CLASSES, f"Class {c} must map to valid canonical 8 class")

    def test_canonical_8_primary_mappings(self):
        expected_mappings = {
            "fev1": "lung_function_fev1",
            "sweat_chloride": "sweat_chloride",
            "exacerbations": "pulmonary_exacerbations",
            "cfqr": "quality_of_life",
            "pro_named_other": "quality_of_life",
            "pro_unnamed": "quality_of_life",
            "microbiology_sputum": "microbiology_cfu",
            "nutrition_growth": "nutritional_bmi",
            "safety_tolerability": "safety_adverse_events",
            "spirometry_other": "other_or_not_stated",
            "lci": "other_or_not_stated",
            "other": "other_or_not_stated",
            "not_stated": "other_or_not_stated",
        }
        for granular, expected_canon in expected_mappings.items():
            self.assertEqual(lexicon.to_canonical_8(granular), expected_canon)

    def test_canonical_classes_cli(self):
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(lexicon.main(["canonical-classes"]), 0)
        output = out.getvalue()
        for canon in lexicon.CANONICAL_8_CLASSES:
            self.assertIn(canon, output)

    def test_tag_entry_populates_canonical_8_classes_and_spans(self):
        entry = {
            "entry_id": "NCT00744159:P1",
            "nct_id": "NCT00744159",
            "measure": "Absolute change from baseline in percent predicted FEV1 through Week 24",
            "description": "",
            "time_frame": "Day 1 to Week 24",
        }
        res = LEX.tag_entry(entry)
        self.assertIn("canonical_8_classes", res)
        self.assertEqual(res["canonical_8_classes"], ["lung_function_fev1"])
        self.assertIn("canonical_8_spans", res)
        self.assertEqual(len(res["canonical_8_spans"]), 1)
        span_info = res["canonical_8_spans"][0]
        self.assertEqual(span_info["canonical_class"], "lung_function_fev1")
        self.assertIsNotNone(span_info["span"])
        start, end = span_info["span"]
        self.assertEqual(entry["measure"][start:end], span_info["text"])

    def test_verbatim_span_integrity_across_synthetic_cases(self):
        for cls, positives, _ in CASES:
            for text in positives:
                entry = {
                    "entry_id": f"SYN:{cls}",
                    "nct_id": "NCT00000000",
                    "measure": text,
                    "description": "",
                    "time_frame": "Week 12",
                }
                res = LEX.tag_entry(entry)
                for s in res["canonical_8_spans"]:
                    if s["field"] and s["span"]:
                        start, end = s["span"]
                        raw_text = entry[s["field"]]
                        self.assertEqual(raw_text[start:end], s["text"])

    def test_span_integrity_failure_on_corrupted_span(self):
        bad_tags = [{
            "class": "fev1",
            "rule_id": "TEST",
            "field": "measure",
            "start": 0,
            "end": 5,
            "text": "DIFFERENT",
        }]
        entry = {"measure": "Exact string"}
        with self.assertRaises(ValueError):
            lexicon.canonical_8_spans(entry, bad_tags)


if __name__ == "__main__":
    unittest.main()
