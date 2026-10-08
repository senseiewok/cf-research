"""Tests for variant_profile.py against the recorded CFTR R31L fixtures. No test reaches the network: every test runs with urlopen replaced by a
guard that fails, and the fetch layer (_http_get) answers from fixtures/. Each rule has a negative control that feeds a bad case and shows the
check catches it.
Usage (from the repo root): python -m unittest discover -s tools/variant_profile -v"""
import copy
import datetime
import io
import json
import sys
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import variant_profile as vp  # noqa: E402

FIX = HERE / "fixtures"
ESEARCH = (FIX / "esearch_CFTR_R31L.json").read_bytes()
ESUMMARY = (FIX / "esummary_54087_53653_35893.json").read_bytes()
VCV = (FIX / "efetch_vcv_54087_is_variationid.xml").read_bytes()
VCV_EMPTY = (FIX / "efetch_vcv_54087.xml").read_bytes()
NOW = datetime.datetime(2026, 10, 7, 12, 0, tzinfo=datetime.timezone.utc)
FOOTER = vp.footer(NOW)


class NetworkUsed(AssertionError):
    pass


def _guard(*a, **k):
    raise NetworkUsed("a test tried to reach the network")


def has_footer(text):
    return "Research tool, not medical advice." in text and "Source: ClinVar via NCBI E-utilities, retrieved 2026-10-07" in text


class Base(unittest.TestCase):
    """Network guarded, pacing disabled, fetch answered from fixtures. self.urls records every URL asked for."""
    vcv = VCV

    def setUp(self):
        vp._last_request[0] = None
        self.urls = []
        for p in (mock.patch("urllib.request.urlopen", _guard), mock.patch.object(vp, "_http_get", self.fake_get),
                  mock.patch.object(vp.time, "sleep", lambda s: None)):
            p.start()
            self.addCleanup(p.stop)

    def fake_get(self, url):
        self.urls.append(url)
        if "esearch.fcgi" in url:
            return ESEARCH
        if "esummary.fcgi" in url:
            return ESUMMARY
        if "efetch.fcgi" in url:
            return self.vcv
        raise AssertionError("unexpected URL " + url)

    def run_main(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        code = vp.main(list(argv), out=out, err=err, now=NOW)
        return code, out.getvalue(), err.getvalue()


def summary_54087():
    return json.loads(ESUMMARY)["result"]["54087"]


class TestParsing(Base):
    def test_profile_fields_from_fixture(self):
        p = vp.build_profile(summary_54087(), VCV)
        self.assertEqual(p["variation_id"], {"value": "54087", "source": "esummary:uid"})
        self.assertEqual(p["germline_classification"]["value"], "Uncertain significance")
        self.assertEqual(p["germline_classification"]["source"], "esummary:germline_classification.description")
        self.assertEqual(p["review_status"]["value"], "reviewed by expert panel")
        self.assertEqual(p["last_evaluated"]["value"], "2017/03/03 00:00")
        self.assertIn("NM_000492.4(CFTR):c.92G>T (p.Arg31Leu)", [h["value"] for h in p["hgvs"]])
        self.assertIn("c.92G>T", [h["value"] for h in p["hgvs"]])
        self.assertEqual(p["number_of_submissions"], {"value": "13", "source": "efetch-vcv:VariationArchive/@NumberOfSubmissions"})
        self.assertEqual(p["supporting_scv_count"]["value"], 13)
        self.assertEqual(len(p["submitters"]), 13)

    def test_submitter_rows_and_counts(self):
        p = vp.build_profile(summary_54087(), VCV)
        cfrow = [r for r in p["submitters"] if r["scv"]["value"] == "SCV001981583"][0]
        self.assertEqual((cfrow["submitter"]["value"], cfrow["classification"]["value"], cfrow["date_last_evaluated"]["value"]),
                         ("CFTR2", "Uncertain significance", "2017-03-03"))
        # 11 + 1 spelled "Uncertain Significance" by one submitter, counted together; 1 "not provided" with no date
        self.assertEqual(p["classification_counts"]["value"], {"Uncertain significance": 12, "not provided": 1})
        staff = [r for r in p["submitters"] if r["classification"]["value"] == "not provided"][0]
        self.assertEqual(staff["date_last_evaluated"]["value"], vp.NOT_STATED)

    def test_every_value_names_its_field(self):
        p = vp.build_profile(summary_54087(), VCV)
        for key, f in p.items():
            items = f if isinstance(f, list) else [f]
            for it in items:
                if key == "submitters":
                    for k in vp.SUBMITTER_FIELDS:
                        self.assertTrue(it[k]["source"].startswith("efetch-vcv:"), k)
                else:
                    self.assertTrue(it["source"], key)

    def test_filter_keeps_only_the_queried_change(self):
        code, out, _ = self.run_main("CFTR", "R31L")
        self.assertEqual(code, 0)
        self.assertIn("54087", out)
        self.assertNotIn("35893", out)  # R31C, also returned by esearch
        self.assertNotIn("53653", out)  # R104fs, also returned by esearch

    def test_three_letter_form_gives_the_same_profile(self):
        code_a, out_a, _ = self.run_main("CFTR", "R31L", "--json")
        vp._last_request[0] = None
        code_b, out_b, _ = self.run_main("CFTR", "Arg31Leu", "--json")
        self.assertEqual(json.loads(out_a)["profiles"], json.loads(out_b)["profiles"])

    def test_json_output(self):
        code, out, _ = self.run_main("CFTR", "R31L", "--json")
        doc = json.loads(out)
        self.assertEqual(code, 0)
        self.assertEqual(doc["profiles"][0]["submitters_disagree"]["value"], False)
        self.assertEqual(doc["catalog_ids"], ["ncbi-eutils", "clinvar"])


class TestNotStated(Base):
    def test_missing_fields_print_not_stated(self):
        s = summary_54087()
        del s["germline_classification"]
        del s["variation_set"]
        p = vp.build_profile(s, None)
        for k in ("germline_classification", "review_status", "last_evaluated"):
            self.assertEqual(p[k]["value"], vp.NOT_STATED, k)
        self.assertEqual(p["hgvs"][0]["value"], vp.NOT_STATED)
        self.assertEqual(p["number_of_submissions"]["value"], vp.NOT_STATED)

    def test_empty_and_placeholder_values_print_not_stated(self):
        s = summary_54087()
        s["germline_classification"]["description"] = ""
        s["germline_classification"]["last_evaluated"] = "1/01/01 00:00"
        p = vp.build_profile(s, None)
        self.assertEqual(p["germline_classification"]["value"], vp.NOT_STATED)
        self.assertEqual(p["last_evaluated"]["value"], vp.NOT_STATED)

    def test_negative_control_present_fields_are_not_replaced(self):
        p = vp.build_profile(summary_54087(), None)
        self.assertNotEqual(p["germline_classification"]["value"], vp.NOT_STATED)

    def test_missing_submitter_attribute_is_not_stated(self):
        body = VCV.replace(b'SubmitterName="CFTR2"', b"")
        rows = vp.parse_vcv(body)[2]
        row = [r for r in rows if r["scv"]["value"] == "SCV001981583"][0]
        self.assertEqual(row["submitter"]["value"], vp.NOT_STATED)

    def test_empty_vcv_set_is_an_api_problem_not_a_profile(self):
        with self.assertRaises(vp.ApiError):
            vp.parse_vcv(VCV_EMPTY)


class TestDisagreement(Base):
    def test_recorded_fixture_does_not_disagree(self):
        # one submitter writes "Uncertain Significance" and one makes no call ("not provided"); neither is a disagreement
        self.assertFalse(vp.build_profile(summary_54087(), VCV)["submitters_disagree"]["value"])

    def test_negative_control_one_differing_classification_raises_the_flag(self):
        body = VCV.replace(b"<GermlineClassification>Uncertain significance</GermlineClassification>",
                           b"<GermlineClassification>Pathogenic</GermlineClassification>", 1)
        self.assertNotEqual(body, VCV)
        p = vp.build_profile(summary_54087(), body)
        self.assertTrue(p["submitters_disagree"]["value"])
        self.assertEqual(p["classification_counts"]["value"]["Pathogenic"], 1)

    def test_flag_in_text_output(self):
        self.vcv = VCV.replace(b"<GermlineClassification>Uncertain significance</GermlineClassification>",
                               b"<GermlineClassification>Likely benign</GermlineClassification>", 1)
        _, out, _ = self.run_main("CFTR", "R31L")
        self.assertIn("submitters disagree", out)
        self.vcv = VCV
        vp._last_request[0] = None
        _, out, _ = self.run_main("CFTR", "R31L")
        self.assertNotIn("submitters disagree", out)
        self.assertIn("submitters do not disagree", out)

    def test_tally_rules(self):
        row = lambda c: {"classification": {"value": c}}  # noqa: E731
        self.assertFalse(vp.tally([row("Pathogenic"), row("pathogenic")])[1])
        self.assertFalse(vp.tally([row("Pathogenic"), row("not provided"), row(vp.NOT_STATED)])[1])
        self.assertTrue(vp.tally([row("Pathogenic"), row("Likely pathogenic")])[1])


class TestInputValidation(Base):
    GOOD = [("CFTR", "R31L"), ("CFTR", "F508del"), ("CFTR", "G542X"), ("CFTR", "R553*"), ("CFTR", "Arg31Leu"), ("CFTR", "Gly542Ter"),
            ("CFTR", "S1255fs"), ("HLA-B", "A2V")]
    BAD = [("John", "Smith"), ("CFTR", "1985-03-02"), ("CFTR", "MRN12345"), ("cftr", "R31L"), ("CFTR", "p.R31L"), ("CFTR", "R31"),
           ("CFTR", "R31L; rm"), ("CFTR", "DOB 02/03/1985"), ("JANE DOE", "R31L"), ("CFTR", "r31l"), ("CFTR", "Xyz31Leu"),
           ("CFTRVERYLONGNAME", "R31L"), ("", "R31L")]

    def test_good_inputs_pass(self):
        for g, c in self.GOOD:
            vp.validate(g, c)

    def test_negative_control_patient_like_inputs_refused_with_exit_2_and_no_request(self):
        for g, c in self.BAD:
            with self.subTest(gene=g, change=c):
                with self.assertRaises(vp.Refused):
                    vp.validate(g, c)
                code, out, err = self.run_main(g, c)
                self.assertEqual(code, 2)
                self.assertIn("refused", err)
                self.assertTrue(has_footer(out))
        self.assertEqual(self.urls, [])

    def test_three_letter_converts_for_the_search(self):
        self.assertEqual(vp.validate("CFTR", "Gly542Ter")[1], "G542X")
        self.assertEqual(vp.validate("CFTR", "R553*")[1], "R553X")


class TestConduct(Base):
    def test_rate_limit_constants(self):
        self.assertLessEqual(vp.MAX_RPS, 1.0)
        self.assertGreaterEqual(vp.MIN_INTERVAL_S, 1.0)
        self.assertLess(vp.MAX_RPS, vp.NCBI_PUBLISHED_RPS)

    def test_second_request_waits(self):
        sleeps, t = [], [100.0]
        vp.fetch("esearch.fcgi", {"db": "clinvar"}, clock=lambda: t[0], sleep=sleeps.append)
        vp.fetch("esearch.fcgi", {"db": "clinvar"}, clock=lambda: t[0], sleep=sleeps.append)
        self.assertEqual(sleeps, [vp.MIN_INTERVAL_S])

    def test_negative_control_no_wait_after_the_interval_passed(self):
        sleeps, t = [], [100.0]
        vp.fetch("esearch.fcgi", {"db": "clinvar"}, clock=lambda: t[0], sleep=sleeps.append)
        t[0] += 2.0
        vp.fetch("esearch.fcgi", {"db": "clinvar"}, clock=lambda: t[0], sleep=sleeps.append)
        self.assertEqual(sleeps, [])

    def test_requests_carry_tool_and_no_email(self):
        self.run_main("CFTR", "R31L")
        self.assertEqual(len(self.urls), 3)
        for u in self.urls:
            self.assertTrue(u.startswith(vp.BASE_URL))
            self.assertIn("tool=" + vp.TOOL_NAME, u)
            self.assertNotIn("email", u)
            self.assertNotIn("api_key", u)
            self.assertNotIn("%40", u)
        self.assertNotIn("@", vp.USER_AGENT)
        self.assertIn("is_variationid=true", self.urls[2])

    def test_negative_control_email_param_refused(self):
        with self.assertRaises(ValueError):
            vp.fetch("esearch.fcgi", {"db": "clinvar", "email": "x"})
        self.assertEqual(self.urls, [])

    def test_no_submitters_makes_two_requests(self):
        code, out, _ = self.run_main("CFTR", "R31L", "--no-submitters")
        self.assertEqual(code, 0)
        self.assertEqual(len(self.urls), 2)
        self.assertIn("Per-submitter rows: not requested", out)

    def test_api_failure_is_not_retried(self):
        calls = []

        def failing(url):
            calls.append(url)
            raise vp.ApiError("boom")
        with mock.patch.object(vp, "_http_get", failing):
            code, out, err = self.run_main("CFTR", "R31L")
        self.assertEqual((code, len(calls)), (3, 1))
        self.assertIn("not retried", err)
        self.assertTrue(has_footer(out))


class TestFooter(Base):
    def test_footer_on_every_outcome(self):
        cases = [("CFTR", "R31L"), ("CFTR", "R31L", "--json"), ("CFTR", "R31L", "--no-submitters"), ("John", "Smith"), ("John", "Smith", "--json"),
                 ("CFTR", "G542X")]  # the last finds no matching record in the R31L fixture: exit 1
        for argv in cases:
            vp._last_request[0] = None
            code, out, _ = self.run_main(*argv)
            with self.subTest(argv=argv):
                self.assertTrue(has_footer(out), out[-300:])
        self.assertEqual(self.run_main("CFTR", "G542X")[0], 1)

    def test_footer_text_exact(self):
        self.assertEqual(FOOTER, "Research tool, not medical advice. A database entry is not an interpretation of any person's genotype; that belongs with "
                                 "their clinical genetics team. Source: ClinVar via NCBI E-utilities, retrieved 2026-10-07 (UTC).")

    def test_negative_control_footer_check_catches_absence(self):
        _, out, _ = self.run_main("CFTR", "R31L")
        self.assertFalse(has_footer(out.replace(FOOTER, "")))


class TestNoNetwork(unittest.TestCase):
    def test_tests_patch_the_fetch(self):
        t = Base("setUp")
        t.setUp()
        try:
            self.assertIs(vp._http_get.__func__ if hasattr(vp._http_get, "__func__") else None, Base.fake_get)
            with self.assertRaises(NetworkUsed):
                __import__("urllib.request").request.urlopen("https://eutils.ncbi.nlm.nih.gov/")
        finally:
            t.doCleanups()

    def test_negative_control_real_fetch_hits_the_guard(self):
        # the unpatched _http_get, run under the guard, is caught before any packet leaves
        with mock.patch("urllib.request.urlopen", _guard):
            with self.assertRaises(NetworkUsed):
                vp._http_get(vp.BASE_URL + "esearch.fcgi?db=clinvar")


if __name__ == "__main__":
    unittest.main()
