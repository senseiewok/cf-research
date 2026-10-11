"""Check leads/leads.json, the lab's public register of research leads, against its rules and sources/catalog.yaml.

A lead is an idea, question or source worth following that is not yet a task. leads/README.md says what each field holds;
this script enforces it:

- the shape: required fields, no unknown fields, the enums (status, source type, tier, effort, guard, needs_person),
  unique ids of the form L-NNNN, ISO dates not after today, updated not before created;
- sources: at least one; a T1 or higher source needs read_on; a conference-private source has ref "private" and no note,
  so the public file carries no detail of it;
- data_needed: every catalog_ids entry exists in sources/catalog.yaml and the stated access equals the catalog's value
  (a single value when the ids share it, or a mapping of id to value); with no catalog ids, access is "not in catalog";
- a done or in-progress lead links a board row or a pull request;
- guard and must_not_claim are non-empty; visibility is public;
- no email-like text, no Windows user path, no NCT trial id (allowlist below, empty), and no quotation longer than
  25 words in any text field.

Output: one `ERROR <id>: ...` line per problem, then `N leads, E errors`. Exit 0 when there are no errors, 1 when any
error, 2 on a usage or input problem (missing file, invalid JSON, wrong top-level shape). Standard library only; no
network. The catalog is read with a small line parser (an entry's `id:` and `access:` lines), not a YAML library. Run
from anywhere: the default paths come from this file's location.
"""
import argparse
import datetime
import json
import re
import sys
from pathlib import Path

STATUSES = ("idea", "scoping", "data-available", "in-progress", "done", "dropped")
SOURCE_TYPES = frozenset({"registry", "literature", "data", "council", "conference-private", "run"})
TIERS = ("T0", "T1", "T2", "T3")
EFFORTS = frozenset({"S", "M", "L"})
GUARDS = frozenset({"computed", "registry-statement-verbatim", "open-access-quote", "lab-position"})
PEOPLE = frozenset({"clinician-reader", "registry-scientist", "maintainer-approval", "none"})
CREATOR_KINDS = frozenset({"human", "model"})
NOT_IN_CATALOG = "not in catalog"
LINKED_STATUSES = frozenset({"done", "in-progress"})
NCT_ALLOWLIST = frozenset()
MAX_QUOTE_WORDS = 25

FIELDS = (
    "id", "title", "question", "why_it_matters", "status", "sources", "data_needed", "effort", "guard",
    "must_not_claim", "needs_person", "links", "created", "updated", "created_by", "visibility",
)
OPTIONAL_FIELDS = ("notes",)
SOURCE_FIELDS = ("type", "ref", "tier", "read_on", "note")
SOURCE_OPTIONAL = ("read_via",)
# How a source was read. A page read through a fetch tool that summarises pages returns the quote as the tool gave it,
# not as the raw page has it: such a source is at most T1 and carries TOOL_NOTE in its note. A search snippet or a memory
# is T0. direct = the page read as it is; local-text-extract = a PDF or file turned into text on this machine;
# lab-run = output of the lab's own script.
READ_VIA = frozenset({"direct", "local-text-extract", "page-summarising-tool", "search-snippet", "memory", "lab-run"})
TOOL_NOTE = "read through a page-summarising tool; quote not checked against the raw page"
DATA_FIELDS = ("catalog_ids", "access", "already_have", "fields_missing")
DATA_OPTIONAL = ("not_in_catalog",)
EFFORT_FIELDS = ("size", "estimate")
GUARD_FIELDS = ("ladder", "note")
LINK_FIELDS = ("board", "proposals", "prs")
CREATOR_FIELDS = ("kind", "model_id")
TOP_LEVEL = frozenset({"schema_version", "note", "leads"})

ID_RE = re.compile(r"^L-\d{4}$")
BOARD_RE = re.compile(r"^T-\d{4}$")
PR_RE = re.compile(r"^(cf-lab|cf-research|cf-skills)#\d+$")
PROPOSAL_RE = re.compile(r"^proposals/[0-9]{4}-[0-9]{2}-[0-9]{2}-[a-z0-9-]+\.md$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
EMAIL_RE = re.compile(r"[^\s@<>()\"',;:]+@[^\s@<>()\"',;:]+\.[A-Za-z]{2,}")
# A drive letter then a Users folder then a name (either slash), or a backslash Users folder; a web address with a
# /users/ path is not a user folder.
WINDOWS_USER_RE = re.compile(r"(?i)(?:\b[a-z]:[\\/]+|\\)users[\\/]+[^\\/\s]+")
NCT_RE = re.compile(r"(?i)\bNCT\s*-?\s*\d{8}\b")
# A quotation: text between double quotes (straight or curly), or between single quotes that open after a space or
# bracket and close before a space or punctuation, so apostrophes inside words are not taken as quote marks.
QUOTE_RES = (
    re.compile(r"\"([^\"]+)\""),
    re.compile(r"“([^”]+)”"),
    re.compile(r"‘([^’]+)’(?![A-Za-z])"),
    re.compile(r"(?:^|(?<=[\s(\[]))'(.+?)'(?=[\s.,;:!?)\]]|$)"),
)
CATALOG_ID_RE = re.compile(r"^  - id:\s*([A-Za-z0-9._-]+)\s*$")
CATALOG_ACCESS_RE = re.compile(r"^    access:\s*([A-Za-z_-]+)\s*(?:#.*)?$")


class InputProblem(Exception):
    """A file is missing, is not valid JSON, or has the wrong top-level shape (exit 2)."""


def default_paths():
    root = Path(__file__).resolve().parents[2]
    return root / "leads" / "leads.json", root / "sources" / "catalog.yaml"


def load_leads(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        raise InputProblem(f"leads file not found: {path}")
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise InputProblem(f"invalid JSON in leads file: {e}")
    if not isinstance(data, dict) or not isinstance(data.get("leads"), list):
        raise InputProblem("leads file must be an object with a 'leads' list")
    return data


def load_catalog_access(path):
    """Map catalog id -> access value, from the entry lines of sources/catalog.yaml."""
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except FileNotFoundError:
        raise InputProblem(f"catalog file not found: {path}")
    except UnicodeDecodeError as e:
        raise InputProblem(f"catalog is not UTF-8: {e}")
    access, current = {}, None
    for line in lines:
        m = CATALOG_ID_RE.match(line)
        if m:
            current = m.group(1)
            access.setdefault(current, None)
            continue
        m = CATALOG_ACCESS_RE.match(line)
        if m and current is not None and access.get(current) is None:
            access[current] = m.group(1)
    if not access:
        raise InputProblem("catalog has no entries (no '  - id:' lines)")
    return access


def iter_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from iter_strings(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from iter_strings(item)


def as_date(value):
    if not isinstance(value, str) or not DATE_RE.match(value):
        return None
    try:
        return datetime.date.fromisoformat(value)
    except ValueError:
        return None


def nonempty_str(value):
    return isinstance(value, str) and bool(value.strip())


def nonempty_str_list(value):
    return isinstance(value, list) and bool(value) and all(nonempty_str(v) for v in value)


def unknown_keys(obj, allowed, where, out):
    for key in obj:
        if key not in allowed:
            out.append(f"{where}: unknown field {key!r}")


def check_text(lead, out):
    for text in iter_strings(lead):
        if EMAIL_RE.search(text):
            out.append("email-like text")
            break
    for text in iter_strings(lead):
        if WINDOWS_USER_RE.search(text):
            out.append("a user folder path")
            break
    for text in iter_strings(lead):
        for m in NCT_RE.finditer(text):
            if re.sub(r"[\s-]", "", m.group(0)).upper() not in NCT_ALLOWLIST:
                out.append(f"trial id {m.group(0)!r} is not on the allowlist")
    for text in iter_strings(lead):
        for rx in QUOTE_RES:
            for m in rx.finditer(text):
                n = len(m.group(1).split())
                if n > MAX_QUOTE_WORDS:
                    out.append(f"a quotation of {n} words; at most {MAX_QUOTE_WORDS}")


def check_date(value, key, today, out):
    d = as_date(value)
    if d is None:
        out.append(f"{key} must be a YYYY-MM-DD date")
    elif d > today:
        out.append(f"{key} {value} is after today ({today.isoformat()})")
    return d


def check_sources(sources, today, out):
    if not isinstance(sources, list) or not sources:
        out.append("sources must be a non-empty list")
        return
    for i, s in enumerate(sources, start=1):
        where = f"sources[{i}]"
        if not isinstance(s, dict):
            out.append(f"{where} must be an object")
            continue
        unknown_keys(s, SOURCE_FIELDS + SOURCE_OPTIONAL, where, out)
        for key in SOURCE_FIELDS:
            if key not in s:
                out.append(f"{where}: missing {key!r}")
        stype, tier = s.get("type"), s.get("tier")
        if "type" in s and stype not in SOURCE_TYPES:
            out.append(f"{where}: type {stype!r} is not one of {sorted(SOURCE_TYPES)}")
        if "ref" in s and not nonempty_str(s.get("ref")):
            out.append(f"{where}: ref must be a non-empty string")
        if "tier" in s and tier not in TIERS:
            out.append(f"{where}: tier {tier!r} is not one of T0..T3")
        if "note" in s and not isinstance(s.get("note"), str):
            out.append(f"{where}: note must be a string")
        read_on = s.get("read_on")
        if read_on is not None:
            check_date(read_on, f"{where}.read_on", today, out)
        elif tier in TIERS and tier != "T0":
            out.append(f"{where}: a {tier} source needs read_on")
        via, note = s.get("read_via"), s.get("note") if isinstance(s.get("note"), str) else ""
        if "read_via" in s and via not in READ_VIA:
            out.append(f"{where}: read_via {via!r} is not one of {sorted(READ_VIA)}")
        if via == "page-summarising-tool":
            if TOOL_NOTE not in note:
                out.append(f"{where}: a source read through a page-summarising tool needs the note {TOOL_NOTE!r}")
            if tier in TIERS and TIERS.index(tier) > 1:
                out.append(f"{where}: a source read through a page-summarising tool is at most T1")
        elif TOOL_NOTE in note:
            out.append(f"{where}: a note that says {TOOL_NOTE!r} needs read_via page-summarising-tool")
        if via in ("search-snippet", "memory") and tier in TIERS and tier != "T0":
            out.append(f"{where}: a source from a {via} is T0")
        if stype == "conference-private":
            if s.get("ref") != "private":
                out.append(f"{where}: a conference-private source has ref 'private' and nothing else")
            if s.get("note"):
                out.append(f"{where}: a conference-private source carries no note in the public file")


def check_data_needed(data, catalog, out):
    if not isinstance(data, dict):
        out.append("data_needed must be an object")
        return
    unknown_keys(data, DATA_FIELDS + DATA_OPTIONAL, "data_needed", out)
    for key in DATA_FIELDS:
        if key not in data:
            out.append(f"data_needed: missing {key!r}")
    ids = data.get("catalog_ids")
    if not isinstance(ids, list) or not all(isinstance(i, str) for i in ids):
        out.append("data_needed.catalog_ids must be a list of catalog ids")
        ids = []
    if len(set(ids)) != len(ids):
        out.append("data_needed.catalog_ids has a duplicate")
    for cid in ids:
        if cid not in catalog:
            out.append(f"data_needed: catalog id {cid!r} is not in sources/catalog.yaml")
    access = data.get("access")
    if not ids:
        if "access" in data and access != NOT_IN_CATALOG:
            out.append(f"data_needed.access must be {NOT_IN_CATALOG!r} when catalog_ids is empty")
    elif isinstance(access, str):
        for cid in ids:
            if cid in catalog and catalog[cid] != access:
                out.append(f"data_needed.access {access!r} differs from the catalog's {catalog[cid]!r} for {cid}"
                           " (use a mapping of id to access when the ids differ)")
    elif isinstance(access, dict):
        if set(access) != set(ids):
            out.append("data_needed.access mapping must name exactly the catalog_ids")
        for cid, value in access.items():
            if cid in catalog and catalog[cid] != value:
                out.append(f"data_needed.access {value!r} differs from the catalog's {catalog[cid]!r} for {cid}")
    elif "access" in data:
        out.append("data_needed.access must be a string or a mapping of catalog id to access")
    if "already_have" in data and not isinstance(data["already_have"], bool):
        out.append("data_needed.already_have must be true or false")
    fm = data.get("fields_missing")
    if "fields_missing" in data and (not isinstance(fm, list) or not all(nonempty_str(f) for f in fm)):
        out.append("data_needed.fields_missing must be a list of strings")
    nic = data.get("not_in_catalog")
    if "not_in_catalog" in data and not nonempty_str_list(nic):
        out.append("data_needed.not_in_catalog, when present, must be a non-empty list of strings")


def check_lead(lead, catalog, today):
    out = []
    unknown_keys(lead, FIELDS + OPTIONAL_FIELDS, "lead", out)
    for key in FIELDS:
        if key not in lead:
            out.append(f"missing required field {key!r}")

    if "id" in lead and (not isinstance(lead["id"], str) or not ID_RE.match(lead["id"])):
        out.append("id must look like L-0001")
    for key in ("title", "question", "why_it_matters"):
        if key in lead and not nonempty_str(lead[key]):
            out.append(f"{key} must be a non-empty string")
    if "notes" in lead and not isinstance(lead["notes"], str):
        out.append("notes must be a string")
    status = lead.get("status")
    if "status" in lead and status not in STATUSES:
        out.append(f"status {status!r} is not one of {list(STATUSES)}")

    if "sources" in lead:
        check_sources(lead["sources"], today, out)
    if "data_needed" in lead:
        check_data_needed(lead["data_needed"], catalog, out)

    effort = lead.get("effort")
    if "effort" in lead:
        if not isinstance(effort, dict):
            out.append("effort must be an object {size, estimate}")
        else:
            unknown_keys(effort, EFFORT_FIELDS, "effort", out)
            if effort.get("size") not in EFFORTS:
                out.append("effort.size must be S, M or L")
            if effort.get("estimate") is not True:
                out.append("effort.estimate must be true (an effort is an estimate until measured on the board)")

    guard = lead.get("guard")
    if "guard" in lead:
        if not isinstance(guard, dict):
            out.append("guard must be an object {ladder, note}")
        else:
            unknown_keys(guard, GUARD_FIELDS, "guard", out)
            ladder = guard.get("ladder")
            if not isinstance(ladder, list) or not ladder:
                out.append("guard.ladder must be a non-empty list")
            else:
                for g in ladder:
                    if g not in GUARDS:
                        out.append(f"guard.ladder has unknown value {g!r}")
            if not nonempty_str(guard.get("note")):
                out.append("guard.note must be a non-empty string")

    if "must_not_claim" in lead and not nonempty_str_list(lead["must_not_claim"]):
        out.append("must_not_claim must be a non-empty list of strings")

    people = lead.get("needs_person")
    if "needs_person" in lead:
        if not isinstance(people, list) or not people:
            out.append("needs_person must be a non-empty list")
        else:
            for p in people:
                if p not in PEOPLE:
                    out.append(f"needs_person has unknown value {p!r}")
            if "none" in people and len(people) > 1:
                out.append("needs_person 'none' stands alone")

    links = lead.get("links")
    if "links" in lead:
        if not isinstance(links, dict):
            out.append("links must be an object {board, proposals, prs}")
            links = {}
        else:
            unknown_keys(links, LINK_FIELDS, "links", out)
            for key, rx in (("board", BOARD_RE), ("proposals", PROPOSAL_RE), ("prs", PR_RE)):
                value = links.get(key, [])
                if not isinstance(value, list):
                    out.append(f"links.{key} must be a list")
                    continue
                for v in value:
                    if not isinstance(v, str) or not rx.match(v):
                        out.append(f"links.{key} has {v!r}, not in the form {rx.pattern}")
        if status in LINKED_STATUSES and not (links.get("board") or links.get("prs")):
            out.append(f"a {status} lead must link a board row or a pull request")

    created = check_date(lead.get("created"), "created", today, out) if "created" in lead else None
    updated = check_date(lead.get("updated"), "updated", today, out) if "updated" in lead else None
    if created and updated and updated < created:
        out.append("updated is before created")

    creator = lead.get("created_by")
    if "created_by" in lead:
        if not isinstance(creator, dict):
            out.append("created_by must be an object {kind, model_id}")
        else:
            unknown_keys(creator, CREATOR_FIELDS, "created_by", out)
            if creator.get("kind") not in CREATOR_KINDS:
                out.append("created_by.kind must be human or model")
            if creator.get("kind") == "model" and not nonempty_str(creator.get("model_id")):
                out.append("created_by.model_id is required when kind is model")

    if "visibility" in lead and lead["visibility"] != "public":
        out.append("visibility must be public; a lead that cannot be public belongs in the private folder, not here")

    check_text(lead, out)
    return out


def check_document(data, catalog, today):
    """Return (lines, n_leads, n_errors)."""
    lines, n_err = [], 0
    for key in data:
        if key not in TOP_LEVEL:
            lines.append(f"ERROR top-level: unknown key {key!r}")
            n_err += 1
    if data.get("schema_version") != 1:
        lines.append(f"ERROR top-level: schema_version must be 1, got {data.get('schema_version')!r}")
        n_err += 1
    seen = set()
    for i, lead in enumerate(data["leads"], start=1):
        label = f"lead {i}"
        if not isinstance(lead, dict):
            lines.append(f"ERROR {label}: must be an object")
            n_err += 1
            continue
        errors = check_lead(lead, catalog, today)
        lid = lead.get("id")
        if isinstance(lid, str) and lid:
            label = lid
            if lid in seen:
                errors.insert(0, f"duplicate id (lead {i})")
            seen.add(lid)
        for msg in errors:
            lines.append(f"ERROR {label}: {msg}")
        n_err += len(errors)
    return lines, len(data["leads"]), n_err


def parse_today(text):
    try:
        return datetime.date.fromisoformat(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a YYYY-MM-DD date: {text!r}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Check leads/leads.json against its rules and the source catalog.")
    ap.add_argument("--leads", help="path to leads.json (default: leads/leads.json in this repository)")
    ap.add_argument("--catalog", help="path to catalog.yaml (default: sources/catalog.yaml in this repository)")
    ap.add_argument("--today", type=parse_today, help="date to check against, YYYY-MM-DD (default: today)")
    try:
        args = ap.parse_args(argv)
    except SystemExit as e:
        return 0 if e.code == 0 else 2
    today = args.today or datetime.date.today()
    d_leads, d_catalog = default_paths()
    try:
        data = load_leads(Path(args.leads) if args.leads else d_leads)
        catalog = load_catalog_access(Path(args.catalog) if args.catalog else d_catalog)
    except InputProblem as e:
        print(f"error: {e}")
        return 2
    lines, n, n_err = check_document(data, catalog, today)
    for line in lines:
        print(line)
    print(f"{n} leads, {n_err} errors")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
