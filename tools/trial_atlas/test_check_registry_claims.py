"""Tests for check_registry_claims.py.

Part (a), the coverage check, runs everywhere, CI included: every registry statement the tools emit must be covered by a claim in
registry_claims.json, and a changed word in a statement must fail. Part (b) checks each claim's quote against saved page text with
tools/claims/check_claims.py. The registry's page text is third-party text and is not committed, so in CI part (b) can only run
against tiny SYNTHETIC evidence files written here from the claims' own quotes; the real check needs a person's saved pages (see the
README). No network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_atlas as ca  # noqa: E402
import check_registry_claims as crc  # noqa: E402
import fetch_snapshot as fs  # noqa: E402


def run(argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = crc.main([str(a) for a in argv])
    return code, out.getvalue()


class CoverageTest(unittest.TestCase):
    def test_the_claims_file_covers_every_emitted_statement(self):
        problems = crc.coverage_problems(crc.load_claims(crc.CLAIMS_PATH))
        self.assertEqual(problems, [])
        code, out = run([])
        self.assertEqual(code, 0, out)
        self.assertIn("part (b) skipped", out)

    def test_a_changed_word_in_a_statement_fails(self):
        for name, value in (("NO_WARRANTY", ca.NO_WARRANTY.replace("liability", "responsibility")),
                            ("SPONSOR_RESPONSIBILITY", ca.SPONSOR_RESPONSIBILITY.replace("See the registry's Disclaimer.", "Trust it.")),
                            ("THIRD_PARTY_COPYRIGHT", ca.THIRD_PARTY_COPYRIGHT.replace("Possessions", "Provinces")),
                            ("KEEP_CURRENT", ca.KEEP_CURRENT.replace("daily", "weekly")),
                            ("DISCLAIMER_LAST_UPDATED", "2023-08-04"),
                            ("DISCLAIMER_URL", "https://clinicaltrials.gov/about-site/disclaimers"),
                            ("TERMS_URL", "https://clinicaltrials.gov/about-site/terms")):
            with self.subTest(name=name), mock.patch.object(ca, name, value):
                problems = crc.coverage_problems(crc.load_claims(crc.CLAIMS_PATH))
                self.assertTrue(any(name in p for p in problems), problems)
                self.assertEqual(run([])[0], 1)

    def test_a_changed_terms_date_fails(self):
        with mock.patch.dict(fs.TERMS, {"terms_last_updated": "2023-02-01"}):
            self.assertTrue(crc.coverage_problems(crc.load_claims(crc.CLAIMS_PATH)))

    def test_a_claim_without_a_quote_or_with_a_bad_field_fails(self):
        claims = crc.load_claims(crc.CLAIMS_PATH)
        for broken in ({**claims[0], "quote": []}, {**claims[0], "quote": "  "}, {k: v for k, v in claims[0].items() if k != "quote"},
                       {**claims[0], "source": "somewhere-else"}, {**claims[0], "scoop": "misspelt scope"}):
            with self.subTest(broken=sorted(broken)), tempfile.TemporaryDirectory() as tmp:
                p = Path(tmp) / "claims.json"
                p.write_text(json.dumps([broken] + claims[1:]), encoding="utf-8")
                code, out = run(["--claims", p])
                self.assertEqual(code, 1, out)


class RequiredClaimTest(unittest.TestCase):
    """A statement another tool (the page generator) relies on, which these tools do not emit, must still have its claim."""

    def test_the_nlm_attribution_claim_is_present_and_exact(self):
        claims = crc.load_claims(crc.CLAIMS_PATH)
        nlm = [c for c in claims if c["id"] == "nlm-developed"]
        self.assertEqual(len(nlm), 1)
        self.assertEqual(nlm[0]["claim"], "ClinicalTrials.gov was developed by the U.S. National Institutes of Health through its "
                                          "National Library of Medicine")
        self.assertEqual(nlm[0]["source"], "ctgov-terms")
        self.assertEqual(crc.coverage_problems(claims), [])

    def test_a_missing_or_changed_required_claim_fails(self):
        claims = crc.load_claims(crc.CLAIMS_PATH)
        without = [c for c in claims if c["id"] != "nlm-developed"]
        self.assertTrue(any("nlm-developed" in p for p in crc.coverage_problems(without)))
        changed = [({**c, "claim": c["claim"].replace("developed", "funded")} if c["id"] == "nlm-developed" else c) for c in claims]
        self.assertTrue(any("nlm-developed" in p for p in crc.coverage_problems(changed)))


class EvidenceTest(unittest.TestCase):
    """Part (b) on SYNTHETIC evidence: each source's text is just the claims' own quotes, under the address header the real saved
    files carry. Not registry page text."""

    def write_evidence(self, tmp, alter=None):
        claims = crc.load_claims(crc.CLAIMS_PATH)
        paths = {}
        for source, url in (("ctgov-terms", ca.TERMS_URL), ("ctgov-disclaimer", ca.DISCLAIMER_URL)):
            quotes = [q for c in claims if c["source"] == source for q in ([c["quote"]] if isinstance(c["quote"], str) else c["quote"])]
            text = f"URL: {url}\nSYNTHETIC evidence for a test, made from the claims' quotes.\n" + "\n".join(quotes) + "\n"
            if alter and alter[0] == source:
                text = text.replace(alter[1], alter[2])
            paths[source] = Path(tmp) / f"{source}.txt"
            paths[source].write_text(text, encoding="utf-8")
        return paths

    def test_part_b_passes_on_matching_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_evidence(tmp)
            code, out = run(["--terms", p["ctgov-terms"], "--disclaimer", p["ctgov-disclaimer"]])
        self.assertEqual(code, 0, out)
        self.assertIn("part (b): 16 claims, 0 failures", out)

    def test_part_b_fails_when_a_quote_is_altered(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_evidence(tmp, alter=("ctgov-terms", "is updated daily", "is updated weekly"))
            code, out = run(["--terms", p["ctgov-terms"], "--disclaimer", p["ctgov-disclaimer"]])
        self.assertEqual(code, 1, out)
        self.assertIn("FAIL updated-daily-keep-current R2", out)

    def test_part_b_fails_when_the_saved_address_differs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_evidence(tmp, alter=("ctgov-disclaimer", ca.DISCLAIMER_URL, "https://example.org/elsewhere"))
            code, out = run(["--terms", p["ctgov-terms"], "--disclaimer", p["ctgov-disclaimer"]])
        self.assertEqual(code, 1, out)
        self.assertIn("names another address", out)

    def test_part_b_needs_both_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.write_evidence(tmp)
            code, out = run(["--terms", p["ctgov-terms"]])
        self.assertEqual(code, 2, out)
        code, out = run(["--terms", Path(tmp) / "gone.txt", "--disclaimer", Path(tmp) / "gone2.txt"])
        self.assertEqual(code, 2, out)


if __name__ == "__main__":
    unittest.main()
