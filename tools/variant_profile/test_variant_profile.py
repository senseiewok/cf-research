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


# --- follow-up fixtures, recorded 2026-10-08: the six variants a 10-variant audit could not score

ES_F508 = (FIX / "esearch_varname_CFTR_F508del.json").read_bytes()
ES_SIX = (FIX / "esearch_varname_CFTR_six_variants.json").read_bytes()  # one OR query over six p. names; the tool's filter must split it
ES_FS = (FIX / "esearch_varname_CFTR_R104fs.json").read_bytes()
ESUM15 = (FIX / "esummary_15_ids_trimmed.json").read_bytes()
VCV_D1152H = (FIX / "efetch_vcv_35867_trimmed.xml").read_bytes()
ES_EMPTY = b'{"esearchresult":{"count":"0","retmax":"0","retstart":"0","idlist":[]}}'


def es_with_count(body, count):
    d = json.loads(body)
    d["esearchresult"]["count"] = str(count)
    return json.dumps(d).encode()


class Audit(Base):
    """Answers from the follow-up fixtures. self.es is a list of esearch answers, used in order."""

    def setUp(self):
        super().setUp()
        self.es = [ES_SIX]
        self.esum = ESUM15
        self.vcv = VCV_D1152H

    def fake_get(self, url):
        self.urls.append(url)
        if "esearch.fcgi" in url:
            return self.es.pop(0)
        if "esummary.fcgi" in url:
            return self.esum
        if "efetch.fcgi" in url:
            return self.vcv
        raise AssertionError("unexpected URL " + url)

    def profile_ids(self, *argv):
        code, out, _ = self.run_main(*argv, "--json")
        doc = json.loads(out)
        return code, doc, [p["variation_id"]["value"] for p in doc["profiles"]]


class TestStopCodons(Audit):
    def test_w1282x_matches_star_in_protein_change(self):
        code, doc, ids = self.profile_ids("CFTR", "W1282X", "--no-submitters")
        self.assertEqual((code, sorted(ids)), (0, ["1300168", "7129", "983867"]))

    def test_y1092x_matches_star_in_protein_change(self):
        code, doc, ids = self.profile_ids("CFTR", "Tyr1092Ter", "--no-submitters")
        self.assertEqual((code, sorted(ids)), (0, ["375475", "38728", "7211"]))

    def test_negative_control_star_and_x_still_need_the_same_change(self):
        self.assertTrue(vp._same_change("W1282*", "W1282X"))
        self.assertFalse(vp._same_change("W1282*", "W1283X"))
        self.assertFalse(vp._same_change("W1282*", "W1282R"))

    def test_three_letter_names(self):
        self.assertEqual([vp.three_letter(c) for c in ("R31L", "G542X", "F508del", "R104fs")],
                         ["p.Arg31Leu", "p.Gly542Ter", "p.Phe508del", "p.Arg104fs"])


class TestSearchTerm(Audit):
    def test_first_term_is_the_hgvs_protein_name(self):
        self.run_main("CFTR", "W1282X", "--no-submitters")
        self.assertIn(vp.urllib.parse.quote_plus('CFTR[gene] AND "p.Trp1282Ter"[varname]'), self.urls[0])
        self.assertEqual(sum("esearch.fcgi" in u for u in self.urls), 1)

    def test_fallback_term_only_when_the_first_finds_nothing(self):
        self.es = [ES_EMPTY, ES_SIX]
        code, doc, ids = self.profile_ids("CFTR", "N1303K", "--no-submitters")
        searches = [u for u in self.urls if "esearch.fcgi" in u]
        self.assertEqual(len(searches), 2)
        self.assertIn(vp.urllib.parse.quote_plus("CFTR[gene] AND N1303K"), searches[1])
        self.assertEqual(doc["query"]["searched"], "CFTR[gene] AND N1303K")

    def test_negative_control_nothing_found_by_either_term_exits_1_after_two_searches(self):
        self.es = [ES_EMPTY, ES_EMPTY]
        code, out, _ = self.run_main("CFTR", "N1303K")
        self.assertEqual(code, 1)
        self.assertEqual(len(self.urls), 2)
        self.assertTrue(has_footer(out))

    def test_fs_variant_found_by_its_p_name(self):
        self.es, self.esum = [ES_FS], ESUMMARY
        code, doc, ids = self.profile_ids("CFTR", "R104fs", "--no-submitters")
        self.assertEqual((code, ids, doc["query"]["match_status"]), (0, ["53653"], "single"))


class TestAmbiguity(Audit):
    def test_n1303k_two_records_both_printed_and_flagged(self):
        code, out, _ = self.run_main("CFTR", "N1303K", "--no-submitters")
        self.assertEqual(code, 0)
        self.assertIn("Match status: ambiguous (2 records match)", out)
        self.assertIn("the tool does not choose", out)
        self.assertIn("4818550", out)
        self.assertIn("7136", out)
        self.assertIn("Match 2 of 2", out)

    def test_r117h_haplotype_record_is_shown_not_dropped(self):
        code, doc, ids = self.profile_ids("CFTR", "R117H", "--no-submitters")
        self.assertEqual(sorted(ids), ["209047", "7109"])
        types = {p["variation_id"]["value"]: p["record_type"]["value"] for p in doc["profiles"]}
        self.assertEqual(types["209047"], "Haplotype")
        self.assertIn("54087", doc["query"]["esearch_ids_not_matched"])  # returned by esearch, not R117H: listed, not profiled

    def test_negative_control_single_match_is_not_flagged(self):
        code, out, _ = self.run_main("CFTR", "R31L", "--no-submitters")
        self.assertIn("Match status: single", out)
        self.assertNotIn("does not choose", out)

    def test_f508del_main_record_matched_by_title(self):
        self.es = [ES_F508]
        code, doc, ids = self.profile_ids("CFTR", "F508del", "--no-submitters")
        # 7105 (the 101-submission record) has an empty protein_change; 634837 is a haplotype (F508del with I1027T) found by its component name
        self.assertEqual(sorted(ids), ["4072070", "634837", "7105"])
        on = {p["variation_id"]["value"]: p["matched_on"] for p in doc["profiles"]}
        self.assertEqual(on, {"4072070": "esummary:protein_change", "7105": "esummary:title",
                              "634837": "esummary:variation_set[0].variation_name"})
        self.assertEqual(doc["query"]["match_status"], "ambiguous (3 records match)")
        self.assertEqual(doc["query"]["esearch_ids_not_matched"], [])

    def test_negative_control_title_route_needs_the_gene_and_the_exact_p_name(self):
        s = json.loads(ESUM15)["result"]["7105"]
        self.assertTrue(vp.matches(s, "CFTR", "F508del"))
        self.assertFalse(vp.matches(s, "CFTR", "F508C"))
        self.assertFalse(vp.matches(s, "SCNN1B", "F508del"))
        s2 = dict(s, title=s["title"].replace("(p.Phe508del)", "p.Phe508del"), variation_set=[])
        self.assertFalse(vp.matches(s2, "CFTR", "F508del"))

    def test_truncated_search_is_said(self):
        self.es = [es_with_count(ES_SIX, 131)]
        code, out, _ = self.run_main("CFTR", "R117H", "--no-submitters")
        self.assertIn("esearch found 131 records but returned 12", out)

    def test_negative_control_complete_search_has_no_truncation_note(self):
        code, out, _ = self.run_main("CFTR", "R117H", "--no-submitters")
        self.assertNotIn("may be incomplete", out)

    def test_one_efetch_for_all_matches(self):
        self.run_main("CFTR", "N1303K")
        efetch = [u for u in self.urls if "efetch.fcgi" in u]
        self.assertEqual(len(efetch), 1)
        self.assertIn("4818550%2C7136", efetch[0])

    def test_too_many_matches_skip_efetch_and_say_so(self):
        with mock.patch.object(vp, "MAX_SUBMITTER_FETCH", 1):
            code, out, _ = self.run_main("CFTR", "N1303K")
        self.assertEqual(sum("efetch.fcgi" in u for u in self.urls), 0)
        self.assertIn("efetch-vcv not requested: more than 1 matches", out)

    def test_record_missing_from_efetch_is_not_stated(self):
        code, doc, _ = self.profile_ids("CFTR", "N1303K")  # this efetch answer holds only 35867
        for p in doc["profiles"]:
            self.assertEqual(p["number_of_submissions"]["value"], vp.NOT_STATED)
            self.assertEqual(p["number_of_submissions"]["source"], "efetch-vcv answer held no record for this id")


class TestOtherClassificationTerms(Audit):
    def test_drug_response_is_stated_and_labelled(self):
        code, doc, ids = self.profile_ids("CFTR", "D1152H")
        p = doc["profiles"][0]
        self.assertEqual((code, ids), (0, ["35867"]))
        self.assertEqual(p["germline_classification"]["value"], "drug response")
        self.assertTrue(p["classification_scale"]["value"].startswith("other ClinVar term"))
        self.assertEqual(p["number_of_submissions"]["value"], "46")
        self.assertEqual(p["classification_counts"]["value"], {"Pathogenic": 41, "Likely pathogenic": 4, "drug response": 1})
        self.assertTrue(p["submitters_disagree"]["value"])

    def test_negative_control_five_tier_terms_are_labelled_five_tier(self):
        for t in ("Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "benign"):
            self.assertEqual(vp.scale(t)["value"], "five-tier pathogenicity term", t)
        for t in ("drug response", "risk factor", "Conflicting classifications of pathogenicity", "Pathogenic, low penetrance"):
            self.assertNotEqual(vp.scale(t)["value"], "five-tier pathogenicity term", t)
        self.assertEqual(vp.scale("")["value"], vp.NOT_STATED)


if __name__ == "__main__":
    unittest.main()
