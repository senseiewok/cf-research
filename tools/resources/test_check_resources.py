"""Tests for check_resources.py: runs it on small fixtures written to temp directories, always with --today 2026-10-09.
Usage: python -m unittest -v (inside tools/resources). Each failure message shows the checker's output."""
import contextlib
import copy
import datetime
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "check_resources.py"
sys.path.insert(0, str(HERE))
import check_resources  # noqa: E402

TODAY = "2026-10-09"
D = datetime.date

CATALOG = """schema_version: 1
sources:
  - id: help-pages
    kind: resource_pages
    access: fetch
    robots: allow
    robots_checked: 2026-09-01
    landing_page: https://WWW.Example.org/help/
  - id: manual-pages
    kind: resource_pages
    access: manual
    robots: unknown
    robots_checked: 2026-09-01
    landing_page: https://manual.example.org/
  - id: stale-pages
    kind: resource_pages
    access: fetch
    robots: allow
    robots_checked: 2026-01-01
    landing_page: https://stale.example.org/
  - id: some-report
    kind: registry_annual_report
    access: fetch
    robots: allow
    robots_checked: 2026-09-01
    landing_page: https://www.example.org/reports/
"""

BASE_ROW = {
    "id": "cf-helpline",
    "name": "CF Helpline",
    "organisation": "Example CF Trust",
    "countries": ["GB"],
    "url": "https://www.example.org/help/line",
    "audience": ["adult_with_cf", "parent_carer_of_child"],
    "category": ["daily_life"],
    "summary": "A free helpline for people with CF and their families, open on weekdays.",
    "who_for": "stated_countries",
    "format": "online",
    "lang": "en-GB",
    "source_id": "help-pages",
    "basis": "person_read_in_browser",
    "summary_attested_by": "maintainer",
    "summary_attested_on": D(2026, 10, 1),
    "confirmed_by": "maintainer",
    "confirmed_on": D(2026, 10, 1),
    "last_checked": D(2026, 10, 1),
    "status": "active",
}


def row(rehash=True, drop=(), **changes):
    """BASE_ROW with changes applied (a value of None removes the key); content_hash recomputed unless rehash=False."""
    r = copy.deepcopy(BASE_ROW)
    for k in drop:
        r.pop(k, None)
    for k, v in changes.items():
        if v is None:
            r.pop(k, None)
        else:
            r[k] = v
    if rehash and r.get("basis") != "lead_only" and "content_hash" not in changes:
        r["content_hash"] = check_resources.compute_hash(r)
    return r


def doc(*rows, **top):
    d = {"schema_version": 1, "resources": list(rows)}
    d.update(top)
    return d


class Base(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp(prefix="cf-resources "))
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        self.catalog = self.dir / "catalog.yaml"
        self.catalog.write_text(CATALOG, encoding="utf-8")
        self.resources = self.dir / "resources.yaml"

    def run_checker(self, data, extra=(), today=TODAY, raw=None):
        if raw is not None:
            self.resources.write_text(raw, encoding="utf-8")
        elif data is not None:
            self.resources.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True), encoding="utf-8")
        argv = ["--resources", str(self.resources), "--catalog", str(self.catalog), "--today", today, *extra]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = check_resources.main(argv)
        return code, buf.getvalue()

    def assertPasses(self, data, **kw):
        code, out = self.run_checker(data, **kw)
        self.assertEqual(code, 0, out)
        self.assertNotIn("ERROR", out, out)
        return out

    def assertError(self, data, text, **kw):
        code, out = self.run_checker(data, **kw)
        self.assertEqual(code, 1, out)
        self.assertIn(text, out, out)
        self.assertTrue(any(line.startswith("ERROR ") and text in line for line in out.splitlines()), out)
        return out


class PassingFixture(Base):
    def test_good_file_passes(self):
        out = self.assertPasses(doc(row(), row(id="second-row", countries=["europe", "global", "IE"], format=None)))
        self.assertTrue(out.strip().endswith("2 rows, 0 errors, 0 warnings"), out)

    def test_lead_only_without_attestation_passes(self):
        self.assertPasses(doc(row(basis="lead_only", drop=("summary_attested_by", "summary_attested_on",
                                                         "confirmed_by", "confirmed_on"))))

    def test_script_fetch_on_fresh_fetch_entry_passes(self):
        self.assertPasses(doc(row(basis="script_fetch")))

    def test_iso_list_has_249_codes(self):
        self.assertEqual(len(check_resources.ISO_COUNTRIES), 249)


class TopLevel(Base):
    def test_schema_version_2(self):
        self.assertError(doc(row(), schema_version=2), "schema_version must be 1")

    def test_missing_schema_version(self):
        self.assertError({"resources": [row()]}, "schema_version must be 1")

    def test_unknown_top_level_key(self):
        self.assertError(doc(row(), extra="x"), "unknown top-level key 'extra'")

    def test_row_not_mapping(self):
        self.assertError(doc(row(), "just text"), "row must be a mapping")


class Ids(Base):
    def test_duplicate_id(self):
        self.assertError(doc(row(), row()), "duplicate id")

    def test_bad_slug(self):
        for bad in ("CF-Helpline", "cf--helpline", "cf_helpline", "-cf", "cf-"):
            with self.subTest(bad=bad):
                self.assertError(doc(row(id=bad)), "id must be a lowercase slug")

    def test_missing_required_field(self):
        self.assertError(doc(row(name=None)), "missing required field 'name'")

    def test_empty_name_and_organisation(self):
        self.assertError(doc(row(name="  ")), "name must be a non-empty string")
        self.assertError(doc(row(organisation="")), "organisation must be a non-empty string")

    def test_unknown_row_field(self):
        self.assertError(doc(row(phone_line="none")), "unknown field 'phone_line'")


class Countries(Base):
    def test_uk_invalid_with_hint(self):
        out = self.assertError(doc(row(countries=["UK"])), "'UK'")
        self.assertIn("use GB", out)

    def test_gb_ok(self):
        self.assertPasses(doc(row(countries=["GB"])))

    def test_regions_ok(self):
        self.assertPasses(doc(row(countries=["europe"])))
        self.assertPasses(doc(row(countries=["global"])))

    def test_lowercase_code_invalid(self):
        self.assertError(doc(row(countries=["gb"])), "'gb'")

    def test_empty_countries(self):
        self.assertError(doc(row(countries=[])), "countries must be a non-empty list")


class Urls(Base):
    def test_http_url(self):
        self.assertError(doc(row(url="http://www.example.org/help/line")), "url must start with https://")

    def test_userinfo(self):
        self.assertError(doc(row(url="https://me@www.example.org/help")), "user info")

    def test_whitespace(self):
        self.assertError(doc(row(url="https://www.example.org/help line")), "whitespace")

    def test_no_host(self):
        self.assertError(doc(row(url="https:///help")), "url has no host")

    def test_host_mismatch_with_landing_page(self):
        self.assertError(doc(row(url="https://other.example.org/help")), "does not match the catalog landing_page host")

    def test_host_compare_is_case_insensitive(self):
        self.assertPasses(doc(row(url="https://www.EXAMPLE.org/help/line")))


class Enums(Base):
    def test_bad_audience(self):
        self.assertError(doc(row(audience=["patients"])), "audience has unknown value 'patients'")

    def test_empty_audience(self):
        self.assertError(doc(row(audience=[])), "audience must be a non-empty list")

    def test_bad_category(self):
        self.assertError(doc(row(category=["fun"])), "category has unknown value 'fun'")

    def test_bad_who_for(self):
        self.assertError(doc(row(who_for="everyone")), "who_for has unknown value")

    def test_bad_format(self):
        self.assertError(doc(row(format="phone")), "format has unknown value")

    def test_format_optional(self):
        self.assertPasses(doc(row(format=None)))

    def test_bad_lang(self):
        for bad in ("English", "EN", "en_GB", "e"):
            with self.subTest(bad=bad):
                self.assertError(doc(row(lang=bad)), "is not a simple BCP 47 tag")

    def test_good_lang(self):
        self.assertPasses(doc(row(lang="en")))
        self.assertPasses(doc(row(lang="pt-BR")))

    def test_bad_status(self):
        self.assertError(doc(row(status="paused")), "status has unknown value")

    def test_bad_basis(self):
        self.assertError(doc(row(basis="guess")), "basis has unknown value")

    def test_list_or_mapping_values_are_errors_not_crashes(self):
        # Found by fuzzing the first version: an unhashable value in a one-word field raised TypeError instead of an error line.
        for field in ("who_for", "format", "basis", "status"):
            for bad in (["x"], {"a": 1}):
                with self.subTest(field=field, bad=bad):
                    self.assertError(doc(row(**{field: bad})), f"{field} has unknown value")


class Catalog(Base):
    def test_unknown_source_id(self):
        self.assertError(doc(row(source_id="nope")), "is not in the catalog")

    def test_wrong_catalog_kind(self):
        self.assertError(doc(row(source_id="some-report")), "source_id must name a resource_pages entry")

    def test_script_fetch_on_manual_entry(self):
        r = row(source_id="manual-pages", url="https://manual.example.org/a", basis="script_fetch")
        self.assertError(doc(r), "basis must be person_read_in_browser or lead_only")

    def test_person_read_on_manual_entry_ok(self):
        self.assertPasses(doc(row(source_id="manual-pages", url="https://manual.example.org/a")))

    def test_script_fetch_with_stale_robots_checked(self):
        r = row(source_id="stale-pages", url="https://stale.example.org/a", basis="script_fetch")
        self.assertError(doc(r), "robots_checked within 183 days")

    def test_robots_checked_183_days_is_fresh(self):
        # 2026-09-01 is 38 days before TODAY; 183 days later than that it is still fresh, 184 it is not.
        self.assertPasses(doc(row(basis="script_fetch")), today="2027-03-03")
        self.assertError(doc(row(basis="script_fetch")), "robots_checked within 183 days", today="2027-03-04")


class Attestation(Base):
    def test_lead_only_with_confirmed_on(self):
        r = row(basis="lead_only", drop=("confirmed_by",))
        self.assertError(doc(r), "a lead_only row cannot be confirmed")

    def test_lead_only_with_content_hash(self):
        r = row(basis="lead_only", drop=("confirmed_by", "confirmed_on"))
        r["content_hash"] = check_resources.compute_hash(r)
        self.assertError(doc(r), "a lead_only row cannot be confirmed")

    def test_missing_attestation(self):
        self.assertError(doc(row(drop=("summary_attested_by", "summary_attested_on"))),
                         "summary_attested_by and summary_attested_on are required together")

    def test_attestation_half_present(self):
        self.assertError(doc(row(drop=("summary_attested_on",))), "required together")

    def test_attested_by_script(self):
        for who in ("script", "Model", "SCRIPT"):
            with self.subTest(who=who):
                self.assertError(doc(row(summary_attested_by=who)), "summary_attested_by must name a person")

    def test_confirmed_by_model(self):
        self.assertError(doc(row(confirmed_by="model")), "confirmed_by must name a person")

    def test_missing_confirmation(self):
        self.assertError(doc(row(drop=("confirmed_by",))), "missing confirmed_by")
        self.assertError(doc(row(drop=("confirmed_on",))), "missing confirmed_on")
        self.assertError(doc(row(rehash=False)), "missing content_hash")

    def test_hash_correct(self):
        self.assertPasses(doc(row()))

    def test_hash_mismatch(self):
        r = row()
        r["summary"] = "A free helpline for people with CF, open on weekdays."
        self.assertError(doc(r), "content changed since confirmation (hash mismatch)")

    def test_hash_not_hex(self):
        self.assertError(doc(row(content_hash="ABC")), "content_hash must be 64 lowercase hex characters")

    def test_compute_hash_known_value(self):
        import hashlib
        r = {"name": "N", "url": "https://a.example.org/", "summary": " S ", "who_for": "not_stated"}
        want = hashlib.sha256("N\nhttps://a.example.org/\n S \nnot_stated".encode("utf-8")).hexdigest()
        self.assertEqual(check_resources.compute_hash(r), want)


class Summary(Base):
    def test_25_words_ok(self):
        self.assertPasses(doc(row(summary=" ".join(["word"] * 25) + ".")))

    def test_26_words_error(self):
        self.assertError(doc(row(summary=" ".join(["word"] * 26) + ".")), "summary has 26 words")

    def test_two_sentences(self):
        self.assertError(doc(row(summary="A helpline. Open on weekdays.")), "summary must be one sentence")

    def test_quote_marks(self):
        for s in ('A so-called "helpline" for families.', "A so-called “helpline” for families."):
            with self.subTest(s=s):
                self.assertError(doc(row(summary=s)), "double quote mark")

    def test_digits(self):
        self.assertError(doc(row(summary="Open 5 days a week for families.")), "digits other than a four-digit year")

    def test_four_digit_year_allowed(self):
        self.assertPasses(doc(row(summary="A guide to benefits, updated for 2026 and still current.")))

    def test_non_year_four_digits(self):
        self.assertError(doc(row(summary="A guide with 1234 tips.")), "digits other than a four-digit year")

    def test_url_like_text(self):
        for s in ("Visit http now.", "See www.example for more.", "Run by example.org staff.",
                  "Run by Example.COM staff.", "Run by example.net staff."):
            with self.subTest(s=s):
                self.assertError(doc(row(summary=s)), "URL-like text")

    def test_at_sign(self):
        self.assertError(doc(row(summary="Ask us @ the clinic.")), "summary must not contain @")


class SummaryEs(Base):
    ES = "Una línea de ayuda para familias, abierta entre semana."

    def test_optional_on_an_english_row(self):
        self.assertPasses(doc(row(summary_es=self.ES)))

    def test_spanish_row_needs_it(self):
        for lang in ("es-MX", "es", "ES-mx"):
            with self.subTest(lang=lang):
                self.assertError(doc(row(lang=lang)), "needs summary_es")

    def test_spanish_row_with_it_passes(self):
        self.assertPasses(doc(row(lang="es-MX", countries=["MX"], summary_es=self.ES)))

    def test_three_letter_language_starting_with_es_is_not_spanish(self):
        self.assertPasses(doc(row(lang="est")))

    def test_same_rules_as_summary(self):
        cases = {
            " ".join(["palabra"] * 26) + ".": "summary_es has 26 words",
            "Una línea. Abierta entre semana.": "summary_es must be one sentence",
            "Abierta 5 días a la semana.": "summary_es must not contain digits",
            'Una "línea" de ayuda.': "summary_es must not contain a double quote mark",
            "Visite www.example para más.": "summary_es must not contain URL-like text",
            "": "summary_es must be a non-empty string",
        }
        for text, message in cases.items():
            with self.subTest(text=text):
                self.assertError(doc(row(summary_es=text)), message)

    def test_email_in_it_is_found(self):
        self.assertError(doc(row(summary_es="Escriba a ayuda@example.org para dudas.")), "phone number or email")

    def test_list_or_mapping_is_an_error_not_a_crash(self):
        for bad in (["x"], {"a": 1}, 5):
            with self.subTest(bad=bad):
                self.assertError(doc(row(summary_es=bad)), "summary_es must be a non-empty string")

    def test_hash_covers_it(self):
        r = row(summary_es=self.ES)
        r["summary_es"] = "Otra frase distinta para familias."
        self.assertError(doc(r), "hash mismatch")

    def test_rows_without_it_keep_their_old_hash(self):
        r = row()
        old = check_resources.hashlib.sha256(
            chr(10).join([r["name"], r["url"], r["summary"], r["who_for"]]).encode("utf-8")).hexdigest()
        self.assertEqual(r["content_hash"], old)
        self.assertNotEqual(row(summary_es=self.ES)["content_hash"], old)


class ContactDetails(Base):
    def test_email_in_summary(self):
        self.assertError(doc(row(summary="Write to help@example.org for advice.")), "phone number or email in a data row")

    def test_email_in_name(self):
        self.assertError(doc(row(name="Helpline help@example.co")), "phone number or email in a data row")

    def test_email_in_list_item(self):
        self.assertError(doc(row(audience=["everyone", "x@y.zz"])), "phone number or email in a data row")

    def test_phone_examples(self):
        for phone in ("+1 844 266 7277", "01 496 2433", "1-800-378-2233", "(01) 496.2433", "8442667277"):
            with self.subTest(phone=phone):
                self.assertError(doc(row(organisation=f"Example Trust {phone}")), "phone number or email in a data row")

    def test_year_is_not_a_phone(self):
        self.assertPasses(doc(row(organisation="Example Trust founded 1998", name="CF Helpline 2026")))

    def test_year_range_is_not_a_phone(self):
        self.assertPasses(doc(row(summary="Registry figures for 2020-2024 explained for families.")))

    def test_short_number_is_not_a_phone(self):
        out = self.run_checker(doc(row(organisation="Example Trust 12 345")))[1]
        self.assertNotIn("phone number or email", out, out)


class Dates(Base):
    def test_future_last_checked(self):
        self.assertError(doc(row(last_checked=D(2026, 10, 10))), "last_checked 2026-10-10 is after today")

    def test_future_attested_on(self):
        self.assertError(doc(row(summary_attested_on=D(2026, 12, 1))), "summary_attested_on 2026-12-01 is after today")

    def test_future_confirmed_on(self):
        self.assertError(doc(row(confirmed_on=D(2027, 1, 1))), "confirmed_on 2027-01-01 is after today")

    def test_invalid_date(self):
        self.assertError(doc(row(last_checked="2026-13-01")), "last_checked must be a valid YYYY-MM-DD date")

    def test_quoted_iso_date_ok(self):
        self.assertPasses(doc(row(last_checked="2026-10-01")))

    def test_retired_without_retired_on(self):
        self.assertError(doc(row(status="retired")), "a retired row needs retired_on")

    def test_retired_with_retired_on_ok(self):
        self.assertPasses(doc(row(status="retired", retired_on=D(2026, 10, 2))))

    def test_retired_on_when_active(self):
        self.assertError(doc(row(retired_on=D(2026, 10, 2))), "retired_on is only allowed when status is retired")


class Freshness(Base):
    def check_age(self, days, categories, status="active"):
        last = D(2026, 10, 9) - datetime.timedelta(days=days)
        return self.run_checker(doc(row(last_checked=last, category=categories, status=status,
                                        **({"retired_on": D(2026, 10, 9)} if status == "retired" else {}))))

    def assertState(self, days, categories, state):
        code, out = self.check_age(days, categories)
        self.assertEqual(code, 0, out)
        due = "WARN cf-helpline: due for a re-check" in out
        expired = "WARN cf-helpline: expired, the page should say not re-checked since" in out
        self.assertEqual((due, expired), (state == "due", state == "expired"), out)
        if state:
            self.assertIn("due for a person to re-check:", out)
            self.assertIn("cf-helpline https://www.example.org/help/line last checked", out)
            self.assertIn("1 warnings", out)
        else:
            self.assertNotIn("due for a person to re-check:", out)
            self.assertIn("0 warnings", out)

    def test_general_thresholds(self):
        self.assertState(180, ["daily_life"], None)
        self.assertState(181, ["daily_life"], "due")
        self.assertState(365, ["daily_life"], "due")
        self.assertState(366, ["daily_life"], "expired")

    def test_money_insurance_thresholds(self):
        self.assertState(90, ["daily_life", "money_insurance"], None)
        self.assertState(91, ["daily_life", "money_insurance"], "due")
        self.assertState(180, ["money_insurance"], "due")
        self.assertState(181, ["money_insurance"], "expired")

    def test_medicines_access_thresholds(self):
        self.assertState(91, ["medicines_access"], "due")
        self.assertState(181, ["medicines_access"], "expired")

    def test_expired_message_names_date(self):
        _, out = self.check_age(400, ["daily_life"])
        last = (D(2026, 10, 9) - datetime.timedelta(days=400)).isoformat()
        self.assertIn(f"not re-checked since {last}", out)
        self.assertIn(f"last checked {last}", out)

    def test_retired_rows_skipped(self):
        code, out = self.check_age(800, ["daily_life"], status="retired")
        self.assertEqual(code, 0, out)
        self.assertNotIn("WARN", out)


class Cli(Base):
    def test_print_hashes(self):
        r1, r2 = row(), row(id="other-row", summary="Another summary for families.")
        code, out = self.run_checker(doc(r1, r2), extra=["--print-hashes"])
        self.assertEqual(code, 0, out)
        self.assertEqual(out.splitlines(), [f"cf-helpline {check_resources.compute_hash(r1)}",
                                            f"other-row {check_resources.compute_hash(r2)}"])

    def test_print_hashes_runs_no_checks(self):
        code, out = self.run_checker(doc(row(url="http://bad"), schema_version=9), extra=["--print-hashes"])
        self.assertEqual(code, 0, out)
        self.assertNotIn("ERROR", out)

    def test_exit_codes_subprocess(self):
        def run(data=None, raw=None, extra=()):
            if raw is not None:
                self.resources.write_text(raw, encoding="utf-8")
            elif data is not None:
                self.resources.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
            p = subprocess.run([sys.executable, str(SCRIPT), "--resources", str(self.resources),
                                "--catalog", str(self.catalog), "--today", TODAY, *extra],
                               capture_output=True, text=True, encoding="utf-8", timeout=60, stdin=subprocess.DEVNULL)
            return p.returncode, p.stdout + p.stderr

        code, out = run(doc(row()))
        self.assertEqual(code, 0, out)
        self.assertEqual(out.strip().splitlines()[-1], "1 rows, 0 errors, 0 warnings")
        code, out = run(doc(row(countries=["UK"])))
        self.assertEqual(code, 1, out)
        self.assertEqual(out.strip().splitlines()[-1], "1 rows, 1 errors, 0 warnings")
        code, out = run(raw="resources: [unclosed")
        self.assertEqual(code, 2, out)
        code, out = run(doc(row()), extra=["--today", "not-a-date"])
        self.assertEqual(code, 2, out)

    def test_missing_resources_file(self):
        code, out = self.run_checker(None)
        self.assertEqual(code, 2, out)
        self.assertIn("error:", out)

    def test_missing_catalog_file(self):
        self.catalog.unlink()
        code, out = self.run_checker(doc(row()))
        self.assertEqual(code, 2, out)

    def test_invalid_yaml(self):
        code, out = self.run_checker(None, raw="resources: [unclosed")
        self.assertEqual(code, 2, out)

    def test_wrong_top_level_shape(self):
        for raw in ("- a\n- b\n", "schema_version: 1\n", "schema_version: 1\nresources: {}\n"):
            with self.subTest(raw=raw):
                code, out = self.run_checker(None, raw=raw)
                self.assertEqual(code, 2, out)

    def test_catalog_without_sources(self):
        self.catalog.write_text("schema_version: 1\n", encoding="utf-8")
        code, out = self.run_checker(doc(row()))
        self.assertEqual(code, 2, out)

    def test_warnings_only_exit_0(self):
        code, out = self.run_checker(doc(row(last_checked=D(2025, 1, 1))))
        self.assertEqual(code, 0, out)
        self.assertIn("WARN", out)


REAL_RESOURCES, REAL_CATALOG = check_resources.get_default_paths()


class RealFiles(unittest.TestCase):
    @unittest.skipUnless(REAL_RESOURCES.exists(), "resources/resources.yaml does not exist yet")
    def test_real_repo_files_pass(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = check_resources.main(["--today", TODAY])
        self.assertEqual(code, 0, buf.getvalue())


if __name__ == "__main__":
    unittest.main()
