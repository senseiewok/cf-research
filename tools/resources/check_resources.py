"""Check resources/resources.yaml (the list of support resources for people with CF) against sources/catalog.yaml.

Every row must name a resource_pages entry of the source catalog whose landing page shares the row's host, carry the
fields the README describes, say how its summary was read (basis), and, unless it is only a lead, record who confirmed it
and a content_hash of name, url, summary, who_for and, when the row has one, summary_es. No row may hold an email address or a phone number. Rows that are
due or expired for a re-check are warnings and are listed for a person at the end.

Output: one line per finding (`ERROR <id>: ...` or `WARN <id>: ...`), then a `due for a person to re-check:` block when
any row is due or expired, then `N rows, E errors, W warnings`. --print-hashes prints `<id> <hash>` per row and runs no
other check.

Exit 0 when there are no errors (warnings allowed), 1 when any error, 2 on an input problem (missing file, invalid YAML,
wrong top-level shape). Standard library plus PyYAML; no network. Run from anywhere: the default paths come from this
file's location.
"""
import argparse
import datetime
import hashlib
import re
import sys
from pathlib import Path
from urllib.parse import urlsplit

import yaml

# ISO 3166-1 alpha-2 officially assigned codes (249).
ISO_COUNTRIES = frozenset("""
AD AE AF AG AI AL AM AO AQ AR AS AT AU AW AX AZ
BA BB BD BE BF BG BH BI BJ BL BM BN BO BQ BR BS BT BV BW BY BZ
CA CC CD CF CG CH CI CK CL CM CN CO CR CU CV CW CX CY CZ
DE DJ DK DM DO DZ
EC EE EG EH ER ES ET
FI FJ FK FM FO FR
GA GB GD GE GF GG GH GI GL GM GN GP GQ GR GS GT GU GW GY
HK HM HN HR HT HU
ID IE IL IM IN IO IQ IR IS IT
JE JM JO JP
KE KG KH KI KM KN KP KR KW KY KZ
LA LB LC LI LK LR LS LT LU LV LY
MA MC MD ME MF MG MH MK ML MM MN MO MP MQ MR MS MT MU MV MW MX MY MZ
NA NC NE NF NG NI NL NO NP NR NU NZ
OM
PA PE PF PG PH PK PL PM PN PR PS PT PW PY
QA
RE RO RS RU RW
SA SB SC SD SE SG SH SI SJ SK SL SM SN SO SR SS ST SV SX SY SZ
TC TD TF TG TH TJ TK TL TM TN TO TR TT TV TW TZ
UA UG UM US UY UZ
VA VC VE VG VI VN VU
WF WS
YE YT
ZA ZM ZW
""".split())
REGIONS = frozenset({"europe", "global"})
COUNTRY_HINTS = {"UK": "use GB", "EL": "use GR", "EU": "use the region value europe"}

AUDIENCES = frozenset({
    "adult_with_cf", "teen_young_adult", "parent_carer_of_child", "partner_family_carer_of_adult",
    "care_team", "researcher", "everyone",
})
CATEGORIES = frozenset({
    "urgent_crisis", "newly_diagnosed_basics", "find_care_team", "mental_health", "medicines_access",
    "money_insurance", "daily_life", "growing_up_transition", "relationships_family_planning",
    "transplant_advanced", "grief_bereavement", "community_peer", "trials_research_registries",
    "advocacy_involved",
})
WHO_FOR = frozenset({"stated_anyone", "stated_countries", "stated_region", "not_stated", "not_checked"})
FORMATS = frozenset({"online", "in_person", "both"})
BASES = frozenset({"person_read_in_browser", "script_fetch", "lead_only"})
STATUSES = frozenset({"active", "retired"})
NON_FETCH_ACCESS = frozenset({"manual", "request", "forbidden"})
MACHINE_NAMES = frozenset({"script", "model"})

REQUIRED = (
    "id", "name", "organisation", "countries", "url", "audience", "category", "summary", "who_for",
    "lang", "source_id", "basis", "last_checked", "status",
)
OPTIONAL = (
    "format", "summary_es", "summary_attested_by", "summary_attested_on", "confirmed_by", "confirmed_on",
    "content_hash", "retired_on",
)
KNOWN_FIELDS = frozenset(REQUIRED + OPTIONAL)
TOP_LEVEL_KEYS = frozenset({"schema_version", "resources"})
CONFIRM_FIELDS = ("confirmed_by", "confirmed_on", "content_hash")

SLUG_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
LANG_RE = re.compile(r"^[a-z]{2,3}(-[A-Za-z0-9]{2,8})*$")
HASH_RE = re.compile(r"^[0-9a-f]{64}$")
EMAIL_RE = re.compile(r"[^\s@<>()\"',;:]+@[^\s@<>()\"',;:]+\.[A-Za-z]{2,}")
# A run of digits joined by up to three spaces, hyphens, dots or parentheses, with an optional leading plus;
# it counts as a phone number when it holds 7 or more digits.
PHONE_RUN_RE = re.compile(r"\+?\(?\d(?:[ \t\-.()]{0,3}\d)+")
ISO_DATE_RUN_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
YEAR_ONLY_RE = re.compile(r"^(?:19|20)\d\d$")
YEAR_RE = re.compile(r"(?<![0-9A-Za-z])(?:19|20)\d\d(?![0-9A-Za-z])")
URLISH_RE = re.compile(r"http|www\.|\.org|\.com|\.net", re.IGNORECASE)
SPANISH_LANG_RE = re.compile(r"^es(-|$)", re.IGNORECASE)
SENTENCE_BREAK_RE = re.compile(r"[.!?][\"')\]]*\s+\S")
QUOTES = ('"', "“", "”")

ROBOTS_MAX_DAYS = 183
DUE_DAYS, EXPIRED_DAYS = 180, 365
FAST_DUE_DAYS, FAST_EXPIRED_DAYS = 90, 180
FAST_CATEGORIES = frozenset({"money_insurance", "medicines_access"})


class InputProblem(Exception):
    """A file is missing, is not valid YAML, or has the wrong top-level shape (exit 2)."""


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Check resources/resources.yaml against the source catalog.")
    parser.add_argument("--resources", default=None, help="Path to resources YAML file")
    parser.add_argument("--catalog", default=None, help="Path to source catalog YAML file")
    parser.add_argument("--today", default=None, type=parse_today, help="Date to check against, YYYY-MM-DD (default: today)")
    parser.add_argument("--print-hashes", action="store_true", help="Print '<id> <hash>' per row and exit 0")
    return parser.parse_args(argv)


def parse_today(text):
    try:
        return datetime.date.fromisoformat(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a YYYY-MM-DD date: {text!r}")


def get_default_paths():
    root = Path(__file__).resolve().parents[2]
    return root / "resources" / "resources.yaml", root / "sources" / "catalog.yaml"


def load_yaml(path, what):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        raise InputProblem(f"{what} file not found: {path}")
    except (yaml.YAMLError, UnicodeDecodeError) as e:
        raise InputProblem(f"Invalid YAML in {what}: {e}")


def load_resources(path):
    data = load_yaml(path, "resources")
    if not isinstance(data, dict):
        raise InputProblem("resources file must be a mapping with a 'resources' list")
    if not isinstance(data.get("resources"), list):
        raise InputProblem("resources file has no 'resources' list")
    return data


def load_catalog(path):
    data = load_yaml(path, "catalog")
    if not isinstance(data, dict) or not isinstance(data.get("sources"), list):
        raise InputProblem("catalog must be a mapping with a 'sources' list")
    entries = {}
    for item in data["sources"]:
        if isinstance(item, dict) and isinstance(item.get("id"), str):
            entries.setdefault(item["id"], item)
    return entries


def compute_hash(row):
    """Lowercase hex SHA-256 of name + newline + url + newline + summary + newline + who_for, values used exactly.

    A row with a summary_es adds a newline and that text at the end, so rows without one keep the hash they had.
    """
    parts = []
    for key in ("name", "url", "summary", "who_for"):
        value = row.get(key, "")
        parts.append(value if isinstance(value, str) else str(value))
    if "summary_es" in row:
        value = row["summary_es"]
        parts.append(value if isinstance(value, str) else str(value))
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()


def is_one_of(value, allowed):
    """True when value is a string in allowed; a list or mapping from a bad data file is never a member (and never a crash)."""
    return isinstance(value, str) and value in allowed


def as_date(value):
    """A date from a YAML date or an ISO date string; None when it is neither."""
    if isinstance(value, datetime.datetime):
        return None
    if isinstance(value, datetime.date):
        return value
    if isinstance(value, str):
        try:
            return datetime.date.fromisoformat(value)
        except ValueError:
            return None
    return None


def iter_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from iter_strings(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from iter_strings(item)


def has_contact_detail(text):
    if EMAIL_RE.search(text):
        return True
    for m in PHONE_RUN_RE.finditer(text):
        run = m.group(0)
        if sum(c.isdigit() for c in run) < 7:
            continue
        if ISO_DATE_RUN_RE.match(run):
            continue  # an ISO date such as 2026-10-09 (a quoted date field) is not a phone number
        groups = re.findall(r"\d+", run)
        if all(YEAR_ONLY_RE.match(g) for g in groups):
            continue  # years only, such as a range 2020-2024
        return True
    return False


def host_of(url):
    try:
        return (urlsplit(url).hostname or "").lower()
    except ValueError:
        return ""


def check_url(url):
    """Return a problem message or None."""
    if not isinstance(url, str) or not url.startswith("https://"):
        return "url must start with https://"
    if any(c.isspace() for c in url):
        return "url must not contain whitespace"
    try:
        parts = urlsplit(url)
        hostname = parts.hostname
        parts.port  # raises ValueError on a malformed port
    except ValueError:
        return "url is malformed"
    if "@" in parts.netloc:
        return "url must not contain user info"
    if not hostname:
        return "url has no host"
    return None


def check_summary(summary, field="summary"):
    """Return a list of problem messages; field names the key in the messages (summary or summary_es)."""
    if not isinstance(summary, str) or not summary.strip():
        return [f"{field} must be a non-empty string"]
    out = []
    words = summary.split()
    if len(words) > 25:
        out.append(f"{field} has {len(words)} words; at most 25")
    if "\n" in summary.strip() or SENTENCE_BREAK_RE.search(summary.strip()):
        out.append(f"{field} must be one sentence")
    if any(q in summary for q in QUOTES):
        out.append(f"{field} must not contain a double quote mark")
    if any(c.isdigit() for c in YEAR_RE.sub("", summary)):
        out.append(f"{field} must not contain digits other than a four-digit year")
    if URLISH_RE.search(summary):
        out.append(f"{field} must not contain URL-like text")
    if "@" in summary:
        out.append(f"{field} must not contain @")
    return out


def check_person(row, key, out):
    value = row.get(key)
    if not isinstance(value, str) or not value.strip():
        out.append(f"{key} must be a non-empty string")
    elif value.strip().lower() in MACHINE_NAMES:
        out.append(f"{key} must name a person, not {value.strip()!r}")


def check_past_date(row, key, today, out):
    value = row.get(key)
    d = as_date(value)
    if d is None:
        out.append(f"{key} must be a valid YYYY-MM-DD date")
    elif d > today:
        out.append(f"{key} {d.isoformat()} is after today ({today.isoformat()})")
    return d


def check_enum_list(row, key, allowed, out):
    value = row.get(key)
    if not isinstance(value, list) or not value:
        out.append(f"{key} must be a non-empty list")
        return []
    for item in value:
        if not isinstance(item, str) or item not in allowed:
            out.append(f"{key} has unknown value {item!r}")
    return [v for v in value if isinstance(v, str)]


def check_row(row, catalog, today):
    """Return (errors, freshness) for one row mapping; freshness is None or (state, last_checked)."""
    errors = []

    for key in row:
        if key not in KNOWN_FIELDS:
            errors.append(f"unknown field {key!r}")
    for key in REQUIRED:
        if key not in row:
            errors.append(f"missing required field {key!r}")

    rid = row.get("id")
    if "id" in row and (not isinstance(rid, str) or not SLUG_RE.match(rid)):
        errors.append("id must be a lowercase slug (letters, digits, single hyphens)")

    for key in ("name", "organisation"):
        if key in row and (not isinstance(row[key], str) or not row[key].strip()):
            errors.append(f"{key} must be a non-empty string")

    if "countries" in row:
        countries = row["countries"]
        if not isinstance(countries, list) or not countries:
            errors.append("countries must be a non-empty list")
        else:
            for c in countries:
                if isinstance(c, str) and (c in ISO_COUNTRIES or c in REGIONS):
                    continue
                hint = COUNTRY_HINTS.get(c) if isinstance(c, str) else None
                if hint is None and isinstance(c, str) and (c.upper() in ISO_COUNTRIES or c.lower() in REGIONS):
                    hint = f"use {c.upper() if c.upper() in ISO_COUNTRIES else c.lower()}"
                msg = f"countries has {c!r}, not an ISO 3166-1 alpha-2 code or europe/global"
                errors.append(msg + (f" ({hint})" if hint else ""))

    url = row.get("url")
    if "url" in row:
        problem = check_url(url)
        if problem:
            errors.append(problem)

    if "audience" in row:
        check_enum_list(row, "audience", AUDIENCES, errors)
    categories = check_enum_list(row, "category", CATEGORIES, errors) if "category" in row else []

    if "summary" in row:
        errors.extend(check_summary(row["summary"]))
    if "summary_es" in row:
        errors.extend(check_summary(row["summary_es"], "summary_es"))

    if "who_for" in row and not is_one_of(row["who_for"], WHO_FOR):
        errors.append(f"who_for has unknown value {row['who_for']!r}")
    if "format" in row and not is_one_of(row["format"], FORMATS):
        errors.append(f"format has unknown value {row['format']!r}")
    if "lang" in row and (not isinstance(row["lang"], str) or not LANG_RE.match(row["lang"])):
        errors.append(f"lang {row['lang']!r} is not a simple BCP 47 tag")
    if isinstance(row.get("lang"), str) and SPANISH_LANG_RE.match(row["lang"]) and "summary_es" not in row:
        errors.append("a Spanish-language row needs summary_es (Mexican Spanish, the lab's own words)")

    entry = None
    if "source_id" in row:
        entry = catalog.get(row["source_id"]) if isinstance(row["source_id"], str) else None
        if entry is None:
            errors.append(f"source_id {row['source_id']!r} is not in the catalog")
        elif entry.get("kind") != "resource_pages":
            errors.append("source_id must name a resource_pages entry")
            entry = None
        else:
            landing = entry.get("landing_page")
            if not isinstance(landing, str) or not host_of(landing):
                errors.append(f"catalog entry {row['source_id']!r} has no landing_page host")
            elif isinstance(url, str) and not check_url(url) and host_of(url) != host_of(landing):
                errors.append(f"url host {host_of(url)} does not match the catalog landing_page host {host_of(landing)}")

    basis = row.get("basis")
    if "basis" in row and not is_one_of(basis, BASES):
        errors.append(f"basis has unknown value {basis!r}")
    if basis == "script_fetch" and entry is not None:
        access = entry.get("access")
        if access in NON_FETCH_ACCESS:
            errors.append(f"the catalog entry's access is {access}: basis must be person_read_in_browser or lead_only")
        elif access != "fetch":
            errors.append(f"script_fetch needs the catalog entry's access to be fetch, not {access!r}")
        else:
            checked = as_date(entry.get("robots_checked"))
            if checked is None or checked > today or (today - checked).days > ROBOTS_MAX_DAYS:
                errors.append(f"script_fetch needs the catalog entry's robots_checked within {ROBOTS_MAX_DAYS} days before today")

    attest_present = [k for k in ("summary_attested_by", "summary_attested_on") if k in row]
    if basis != "lead_only" or attest_present:
        if len(attest_present) < 2:
            errors.append("summary_attested_by and summary_attested_on are required together")
        if "summary_attested_by" in row:
            check_person(row, "summary_attested_by", errors)
        if "summary_attested_on" in row:
            check_past_date(row, "summary_attested_on", today, errors)

    if basis == "lead_only":
        if any(k in row for k in CONFIRM_FIELDS):
            errors.append("a lead_only row cannot be confirmed")
    elif is_one_of(basis, BASES):
        for k in CONFIRM_FIELDS:
            if k not in row:
                errors.append(f"missing {k} (required unless basis is lead_only)")
        if "confirmed_by" in row:
            check_person(row, "confirmed_by", errors)
        if "confirmed_on" in row:
            check_past_date(row, "confirmed_on", today, errors)
        if "content_hash" in row:
            h = row["content_hash"]
            if not isinstance(h, str) or not HASH_RE.match(h):
                errors.append("content_hash must be 64 lowercase hex characters")
            elif h != compute_hash(row):
                errors.append("content changed since confirmation (hash mismatch)")

    last_checked = check_past_date(row, "last_checked", today, errors) if "last_checked" in row else None

    status = row.get("status")
    if "status" in row and not is_one_of(status, STATUSES):
        errors.append(f"status has unknown value {status!r}")
    if status == "retired":
        if "retired_on" not in row:
            errors.append("a retired row needs retired_on")
        else:
            check_past_date(row, "retired_on", today, errors)
    elif "retired_on" in row:
        errors.append("retired_on is only allowed when status is retired")

    for key, value in row.items():
        if key == "content_hash":
            continue  # a hex digest, checked above; its digit runs are not phone numbers
        if any(has_contact_detail(s) for s in iter_strings(value)):
            errors.append("phone number or email in a data row")
            break

    freshness = None
    if status == "active" and last_checked is not None and last_checked <= today:
        fast = any(c in FAST_CATEGORIES for c in categories)
        due_days, expired_days = (FAST_DUE_DAYS, FAST_EXPIRED_DAYS) if fast else (DUE_DAYS, EXPIRED_DAYS)
        age = (today - last_checked).days
        if age > expired_days:
            freshness = ("expired", last_checked)
        elif age > due_days:
            freshness = ("due", last_checked)
    return errors, freshness


def check_document(data, catalog, today):
    """Return (lines, due_lines, n_rows, n_errors, n_warnings)."""
    lines, due_lines = [], []
    n_err = n_warn = 0

    for key in data:
        if key not in TOP_LEVEL_KEYS:
            lines.append(f"ERROR top-level: unknown top-level key {key!r}")
            n_err += 1
    if data.get("schema_version") != 1:
        lines.append(f"ERROR top-level: schema_version must be 1, got {data.get('schema_version')!r}")
        n_err += 1

    rows = data["resources"]
    seen = set()
    for i, row in enumerate(rows, start=1):
        label = f"row {i}"
        if not isinstance(row, dict):
            lines.append(f"ERROR {label}: row must be a mapping")
            n_err += 1
            continue
        rid = row.get("id")
        if isinstance(rid, str) and rid.strip():
            label = rid
        errors, freshness = check_row(row, catalog, today)
        if isinstance(rid, str) and rid.strip():
            if rid in seen:
                errors.insert(0, f"duplicate id (row {i})")
            seen.add(rid)
        for msg in errors:
            lines.append(f"ERROR {label}: {msg}")
        n_err += len(errors)
        if freshness:
            state, checked = freshness
            if state == "expired":
                lines.append(f"WARN {label}: expired, the page should say not re-checked since {checked.isoformat()}")
            else:
                lines.append(f"WARN {label}: due for a re-check")
            n_warn += 1
            due_lines.append(f"{label} {row.get('url')} last checked {checked.isoformat()}")
    return lines, due_lines, len(rows), n_err, n_warn


def main(argv=None):
    args = parse_args(argv)
    today = args.today or datetime.date.today()

    default_resources, default_catalog = get_default_paths()
    resources_path = Path(args.resources) if args.resources else default_resources
    catalog_path = Path(args.catalog) if args.catalog else default_catalog

    try:
        data = load_resources(resources_path)
    except InputProblem as e:
        print(f"error: {e}")
        return 2

    if args.print_hashes:
        for i, row in enumerate(data["resources"], start=1):
            if isinstance(row, dict):
                rid = row.get("id") if isinstance(row.get("id"), str) else f"row-{i}"
                print(f"{rid} {compute_hash(row)}")
        return 0

    try:
        catalog = load_catalog(catalog_path)
    except InputProblem as e:
        print(f"error: {e}")
        return 2

    lines, due_lines, n_rows, n_err, n_warn = check_document(data, catalog, today)
    for line in lines:
        print(line)
    if due_lines:
        print("due for a person to re-check:")
        for line in due_lines:
            print(line)
    print(f"{n_rows} rows, {n_err} errors, {n_warn} warnings")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
