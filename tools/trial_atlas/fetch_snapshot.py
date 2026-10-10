#!/usr/bin/env python3
"""Fetch one dated, hashed snapshot of ClinicalTrials.gov API v2 studies for the trial endpoint atlas (board row T-0130).

Every request goes through the evidence skill's Client (cf-skills, cf-evidence-loop), imported from the sibling checkout exactly as
tools/sources/fetch_sources.py imports it, so its rules apply unchanged: the catalog gate (`clinicaltrials-gov` must be `access: api`
with terms and a rate), 1 request per second, a ceiling of 200 requests per process, no redirect followed, a body cap, one honest
project user agent. This tool adds: a contact is required (it goes into the User-Agent and nowhere else), every page is kept as
received, the registry's dataTimestamp is recorded, and the run fails closed when an answer does not look as expected.

What it writes (into a new folder; it never overwrites one):
    version.json                 the answer of the `version` endpoint (holds dataTimestamp)
    condition/page-0001.json ... the pages of the main route (the condition search)
    term/page-0001.json ...      the pages of the recall route (a term search), unless --route condition
    manifest.json                written last, only when every check passed: the query strings verbatim, the request count,
                                 dataTimestamp, the sha256 of each file and of the whole snapshot. No contact is stored.
A run that stops leaves INCOMPLETE.txt with the reason and no manifest, so nothing downstream reads it.

The default folder is under sources/downloads/, which git ignores: registry pages are not committed by this tool.

NOT VERIFIED against the real registry (nothing here has been run against it). Each item is configuration, with the default
marked "to verify on the first real run" in DEFAULTS below: the largest pageSize, the request parameter that carries the page token,
the query parameters, the `fields` piece names, and whether totalCount comes back on the first page only.

Usage (the contact comes ONLY from the environment variable EVIDENCE_CONTACT, so it never sits on a command line or in shell
history; there is no --contact option, and without the variable the tool does not run):
    python fetch_snapshot.py [--route both|condition|term] [--out DIR] [--page-size N] [--max-pages N]
    python fetch_snapshot.py --dry-run       print the planned requests; no socket is opened
A real run is a network call: it needs a person's approval first (design, build step 2).
Exit 0 written and checked; 1 the run stopped (reason printed, INCOMPLETE.txt left; also on any unexpected error from the client);
2 refused before any request (no contact, a contact the client refuses, evidence skill absent, catalog gate).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from snapshot import VERSION, manifest_digest, sha256_file, valid_nct  # noqa: E402

REPO_ROOT = HERE.parents[1]
CATALOG = REPO_ROOT / "sources" / "catalog.yaml"
DEFAULT_OUT_PARENT = REPO_ROOT / "sources" / "downloads" / "trial_atlas"
TOOL = "trial_atlas/fetch_snapshot 0.1.0"
SOURCE_ID = "clinicaltrials-gov"

# Every value here is an assumption about the registry's behaviour until a first real, approved run confirms it.
DEFAULTS = {
    "studies_path": "studies",               # to verify on the first real run
    "version_path": "version",               # stated on the registry's data refresh page (design); to verify on the first real run
    "page_size": 100,                        # to verify on the first real run: the maximum pageSize is not stated in the documentation
    "page_size_param": "pageSize",           # to verify on the first real run (the evidence skill's provider sends it)
    "page_token_param": "pageToken",         # to verify on the first real run: the request key that carries the next token
    "next_token_field": "nextPageToken",     # to verify on the first real run: the response key that holds it
    "total_field": "totalCount",             # to verify on the first real run (the evidence skill's provider reads it)
    "count_total_param": "countTotal",       # to verify on the first real run: asks for totalCount, sent on the first page only
    "studies_field": "studies",
    "data_timestamp_field": "dataTimestamp", # named on the registry's data refresh page (design)
    # to verify on the first real run: the piece names. The evidence skill's provider uses NCTId, BriefTitle, OverallStatus, Phase,
    # StudyType, StartDate, LeadSponsorName, Condition, PrimaryOutcomeMeasure and HasResults; the others are assumed by analogy.
    "fields": ("NCTId,BriefTitle,OverallStatus,StartDate,StartDateType,StudyFirstPostDate,StudyType,Phase,Condition,Keyword,"
               "PrimaryOutcomeMeasure,PrimaryOutcomeDescription,PrimaryOutcomeTimeFrame,LeadSponsorName,HasResults"),
    "max_pages": 150,                        # our own cap per run, under the evidence client's 200-request ceiling
}

# The two retrieval routes (design, decision 3). The parameter names and the term syntax are to verify on the first real run.
ROUTES = {
    "condition": {"query.cond": "cystic fibrosis"},
    "term": {"query.term": "\"cystic fibrosis\" OR mucoviscidosis"},
}

TERMS = {
    "source": "ClinicalTrials.gov",
    "terms_last_updated": "2023-01-31",
    "on_any_publication": [
        "attribute the source as ClinicalTrials.gov",
        "display the date the data were processed by ClinicalTrials.gov (dataTimestamp)",
        "state every modification made to the content, with a complete description",
        "do not represent the database or any part of it as other than a United States Government database",
        "keep the data current, or state the snapshot date prominently (a frozen snapshot cannot be current at all times)",
        "some content may be subject to third-party copyright; quote only what each count needs",
    ],
    "read_from": "the registry's Terms and Conditions page, as summarised in proposals/2026-10-10-trial-endpoint-atlas.md",
}

_EMAIL_LIKE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+")


class Refused(RuntimeError):
    """A run this tool will not start, with the reason. Nothing was requested."""


class SnapshotStopped(RuntimeError):
    """The run stopped part-way because an answer did not look as expected, or a limit was reached. No manifest is written."""


# ------------------------------------------------------------------ the evidence client, from the sibling checkout

def _evidence_dirs() -> list[Path]:
    env = os.environ.get("EVIDENCE_SKILL_DIR")
    return ([Path(env)] if env else []) + [REPO_ROOT.parent / "cf-skills" / ".claude" / "skills" / "cf-evidence-loop" / "scripts"]


def load_evidence_http():
    """The evidence skill's http module (it holds the Client and every network rule). Its catalog gate reads EVIDENCE_CATALOG; when
    that is not set it is pointed at this repository's sources/catalog.yaml, the single record of what each source permits."""
    os.environ.setdefault("EVIDENCE_CATALOG", str(CATALOG))
    for d in _evidence_dirs():
        if (d / "evidence" / "http.py").is_file():
            sys.path.insert(0, str(d))
            try:
                from evidence import http as module  # noqa: PLC0415
            except ImportError as exc:
                raise Refused(f"the evidence skill could not be imported ({exc}); install tools/sources/requirements.txt") from exc
            finally:
                sys.path.remove(str(d))
            return module
    raise Refused("the evidence skill (cf-skills, cf-evidence-loop) is not available; clone cf-skills next to this repository or set "
                  "EVIDENCE_SKILL_DIR. This tool makes no request without its rules.")


CONTACT_REFUSED = "the contact was refused by the evidence client: not a plain email address"


def resolve_contact() -> str:
    """The contact, from EVIDENCE_CONTACT only (never a command-line argument)."""
    contact = os.environ.get("EVIDENCE_CONTACT") or ""
    if not contact.strip():
        raise Refused("no contact: set the environment variable EVIDENCE_CONTACT. The registry client sends it in the User-Agent; "
                      "this tool does not run without one")
    return contact.strip()


def make_client(contact: str, evidence=None, *, dry_run: bool = False):
    ev = evidence or load_evidence_http()
    if dry_run:
        os.environ["EVIDENCE_DRY_RUN"] = "1"
    try:
        client = ev.Client(contact=contact)
    except ValueError:
        # A fixed message: the client's own message quotes the value, and the value is never echoed.
        raise Refused(CONTACT_REFUSED) from None
    return client, ev


# ------------------------------------------------------------------ one run

def _status(fetch) -> str:
    s = getattr(fetch, "status", None)
    return str(getattr(s, "value", s))


def _request(client, path: str, params: dict | None):
    """One request through the client. The client's gate errors become Refused; its budget error stops the run."""
    try:
        return client.get(SOURCE_ID, path, params=params)
    except (PermissionError, KeyError) as exc:          # catalog.AccessDenied, UnknownSource, PaperworkMissing, HostInCooldown
        raise Refused(f"refused by the evidence client's gate: {type(exc).__name__}: {exc}") from exc
    except RuntimeError as exc:                           # http.BudgetExhausted; any other RuntimeError is not ours to explain
        if type(exc).__name__ != "BudgetExhausted":
            raise
        raise SnapshotStopped(f"request budget exhausted ({type(exc).__name__}: {exc}); start a new process deliberately, "
                              "do not loop") from exc


def _write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def fetch_version(client, cfg: dict) -> dict:
    f = _request(client, cfg["version_path"], None)
    if getattr(f, "dry_run", False):
        print(f"planned  {f.url}")
        return {}
    if _status(f) != "found" or not isinstance(f.data, dict):
        raise SnapshotStopped(f"the version endpoint did not answer as expected (status {_status(f)}, http {getattr(f, 'http_status', None)})")
    if not f.data.get(cfg["data_timestamp_field"]):
        raise SnapshotStopped(f"the version answer has no {cfg['data_timestamp_field']}; the data refresh date cannot be recorded")
    return f.data


def fetch_route(client, name: str, base_params: dict, cfg: dict, out_dir: Path) -> dict:
    """Page through one route. Fails closed on: no studies list, a study without a valid NCT id, no total on the first page, a total
    that changes, an empty page that still names a next token, a repeated token, more pages than max_pages, duplicate studies, and a
    final count that differs from the total (which is also how a missing next token shows)."""
    first_params = {**base_params, cfg["page_size_param"]: cfg["page_size"], cfg["count_total_param"]: "true", "fields": cfg["fields"]}
    route = {"name": name, "query": urlencode(first_params), "params": first_params, "pages": [], "total_count": None,
             "studies_received": 0}
    token, seen_tokens, seen_ids = None, set(), set()
    for n in range(1, cfg["max_pages"] + 1):
        params = dict(first_params)
        if token is not None:
            params.pop(cfg["count_total_param"], None)
            params[cfg["page_token_param"]] = token
        f = _request(client, cfg["studies_path"], params)
        if getattr(f, "dry_run", False):
            print(f"planned  {f.url}")
            return route
        if _status(f) != "found":
            raise SnapshotStopped(f"{name} page {n}: status {_status(f)} (http {getattr(f, 'http_status', None)}); not retried")
        data = f.data
        studies = data.get(cfg["studies_field"]) if isinstance(data, dict) else None
        if not isinstance(studies, list):
            raise SnapshotStopped(f"{name} page {n}: the answer has no '{cfg['studies_field']}' list; the API may have changed")
        for s in studies:
            nct = ((s or {}).get("protocolSection") or {}).get("identificationModule", {}).get("nctId") if isinstance(s, dict) else None
            if not valid_nct(nct):
                raise SnapshotStopped(f"{name} page {n}: a study has no valid protocolSection.identificationModule.nctId")
            if nct in seen_ids:
                raise SnapshotStopped(f"{name} page {n}: {nct} was already received on an earlier page; paging is not stable")
            seen_ids.add(nct)
        total = data.get(cfg["total_field"])
        if n == 1:
            if not isinstance(total, int) or isinstance(total, bool) or total < 0:
                raise SnapshotStopped(f"{name}: the first page has no usable '{cfg['total_field']}' (was {cfg['count_total_param']} "
                                      "honoured?); the received count cannot be checked, so the run stops")
            route["total_count"] = total
        elif total is not None and total != route["total_count"]:
            raise SnapshotStopped(f"{name} page {n}: '{cfg['total_field']}' changed from {route['total_count']} to {total} during "
                                  "the run (a data refresh?); start again after the refresh")
        rel = f"{name}/page-{n:04d}.json"
        _write_json(out_dir / rel, data)
        route["pages"].append({"file": rel, "studies": len(studies), "url": getattr(f, "url", None),
                               "query": urlencode(params)})
        route["studies_received"] += len(studies)
        nxt = data.get(cfg["next_token_field"])
        if nxt in (None, ""):
            break
        if not isinstance(nxt, str):
            raise SnapshotStopped(f"{name} page {n}: '{cfg['next_token_field']}' is not a string")
        if not studies:
            raise SnapshotStopped(f"{name} page {n}: an empty page names a next token; stopping instead of looping")
        if nxt in seen_tokens:
            raise SnapshotStopped(f"{name} page {n}: the next token repeats; stopping instead of looping")
        seen_tokens.add(nxt)
        token = nxt
    else:
        raise SnapshotStopped(f"{name}: more than {cfg['max_pages']} pages; raise --max-pages deliberately or narrow the query")
    if route["studies_received"] != route["total_count"]:
        hint = (f" and the last page named no '{cfg['next_token_field']}'" if route["studies_received"] < route["total_count"] else "")
        raise SnapshotStopped(f"{name}: received {route['studies_received']} studies but '{cfg['total_field']}' is "
                              f"{route['total_count']}{hint}; the snapshot is not complete")
    return route


def build_manifest(out_dir: Path, version: dict, routes: list[dict], cfg: dict, attempts: int, project_ua: str,
                   catalog_entry: dict | None) -> dict:
    files = {}
    for p in sorted(out_dir.rglob("*")):
        rel = p.relative_to(out_dir).as_posix()
        if p.is_file() and rel != "manifest.json":
            files[rel] = sha256_file(p)
    manifest = {
        "tool": TOOL,
        "source_id": SOURCE_ID,
        "base_url": (catalog_entry or {}).get("base_url"),
        "terms_url": (catalog_entry or {}).get("terms_url"),
        "fetched_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_timestamp": version.get(cfg["data_timestamp_field"]),
        "api_version": version.get("apiVersion"),
        "user_agent_product": project_ua,
        "contact_recorded": False,
        "requests": attempts,
        "config": {k: v for k, v in cfg.items()},
        "config_verified": False,
        "routes": [{k: r[k] for k in ("name", "query", "params", "total_count", "studies_received", "pages")} for r in routes],
        "terms": TERMS,
        "modifications": "none at fetch: each page is the parsed JSON answer written back as UTF-8 JSON (the evidence client returns "
                         "parsed data, not bytes), so key order and values are kept and whitespace is not",
        "files": files,
    }
    manifest["snapshot_sha256"] = manifest_digest(manifest)     # covers the whole manifest, not only the file list
    return manifest


def run(client, out_dir: Path, *, routes: list[str], cfg: dict | None = None, project_ua: str = "", contact: str = "",
        catalog_entry: dict | None = None, route_params: dict | None = None) -> dict | None:
    """Fetch into out_dir (which must not exist). Returns the manifest, or None for a dry run. On SnapshotStopped or Refused after
    the folder was made, INCOMPLETE.txt holds the reason and no manifest is written. route_params replaces ROUTES (used to build the
    synthetic non-CF control snapshot, whose query must not claim to be a CF search)."""
    cfg = {**DEFAULTS, **(cfg or {})}
    params_by_route = route_params or ROUTES
    unknown = [r for r in routes if r not in params_by_route]
    if unknown:
        raise Refused(f"unknown route(s): {', '.join(unknown)}")
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise Refused(f"{out_dir.name} already exists; a snapshot is never overwritten")
    dry = bool(getattr(client, "dry_run", False))
    if not dry:
        out_dir.mkdir(parents=True)
    try:
        version = fetch_version(client, cfg)
        done = [fetch_route(client, r, params_by_route[r], cfg, out_dir) for r in routes]
        if dry:
            return None
        _write_json(out_dir / VERSION, version)
        manifest = build_manifest(out_dir, version, done, cfg, client.accounting.attempts, project_ua, catalog_entry)
        text = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
        if contact and contact in text:
            raise SnapshotStopped("the contact would have been written into the manifest; refused")
        if _EMAIL_LIKE.search(text):
            raise SnapshotStopped("something shaped like an email address would have been written into the manifest; refused")
        (out_dir / "manifest.json").write_text(text, encoding="utf-8")
        return manifest
    except Exception as exc:  # noqa: BLE001 - every failure leaves INCOMPLETE.txt, then propagates unchanged
        if out_dir.exists():
            attempts = getattr(getattr(client, "accounting", None), "attempts", "?")
            # Our own errors carry a reason we wrote; anything else is named by its type only, so no text from elsewhere is copied.
            why = str(exc) if isinstance(exc, (SnapshotStopped, Refused)) else f"unexpected error from the client: {type(exc).__name__}"
            (out_dir / "INCOMPLETE.txt").write_text(f"stopped: {why}\nrequests made: {attempts}\nno manifest was written\n",
                                                    encoding="utf-8")
        raise


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--route", choices=("both", "condition", "term"), default="both")
    ap.add_argument("--out", type=Path, help="a new folder for the snapshot (default: under sources/downloads/trial_atlas/)")
    ap.add_argument("--page-size", type=int, default=DEFAULTS["page_size"])
    ap.add_argument("--max-pages", type=int, default=DEFAULTS["max_pages"])
    ap.add_argument("--dry-run", action="store_true", help="print the planned requests; open no socket")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    if a.page_size < 1 or a.max_pages < 1:
        print("error: --page-size and --max-pages must be positive", file=sys.stderr)
        return 2
    try:
        contact = resolve_contact()
        client, ev = make_client(contact, dry_run=a.dry_run)
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    try:
        from evidence import catalog as ev_catalog  # noqa: PLC0415  (loaded with the client above)
        entry = ev_catalog.entry(SOURCE_ID)
    except Exception:  # noqa: BLE001 - the gate inside the client still decides; this is only for the manifest
        entry = None
    routes = ["condition", "term"] if a.route == "both" else [a.route]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = a.out or DEFAULT_OUT_PARENT / f"snapshot-{stamp}"
    try:
        manifest = run(client, out, routes=routes, cfg={"page_size": a.page_size, "max_pages": a.max_pages},
                       project_ua=ev.PROJECT_UA, contact=contact, catalog_entry=entry)
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    except SnapshotStopped as exc:
        print(f"STOPPED: {exc}", file=sys.stderr)
        print(client.accounting.summary(), file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 - INCOMPLETE.txt is already written; name the error by type only
        print(f"STOPPED: unexpected error from the client: {type(exc).__name__}", file=sys.stderr)
        return 1
    print(client.accounting.summary())
    if manifest is None:
        print("dry run: nothing was requested or written")
        return 0
    print(f"snapshot: {out.name}  sha256 {manifest['snapshot_sha256']}")
    print(f"dataTimestamp: {manifest['data_timestamp']}  requests: {manifest['requests']}")
    for r in manifest["routes"]:
        print(f"route {r['name']}: {r['studies_received']} studies received, total {r['total_count']}, {len(r['pages'])} page(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
