"""Read a trial-atlas snapshot written by fetch_snapshot.py, check it against its manifest, and turn its studies into plain records.

A snapshot is a folder:  manifest.json, version.json and one folder per retrieval route (condition/, term/) of page-NNNN.json files.
The manifest holds the sha256 of every other file, and snapshot_sha256: the sha256 of a canonical form of the whole manifest without
that one key (manifest_digest), so the routes, the query strings, the totals, dataTimestamp and the file list are all covered.
Loading refuses a snapshot unless:
  - snapshot_sha256 matches the manifest as it is on disk;
  - the manifest's file list equals the files under the folder exactly, as spelled on disk (no alias such as './', another letter
    case or a trailing dot), and every listed file matches its hash;
  - every page a route names, and version.json, is a listed file, given as a relative path inside the folder (no absolute path,
    no '..', nothing that resolves outside the folder); no page is listed twice and no route name twice;
  - every study has an NCT id of the form NCT followed by 8 digits, no study appears twice in a route, each page holds the number
    of studies the manifest says, and each route's studies add up to its studies_received and its total_count.
So a count is only ever computed from bytes the manifest names. The hash is a consistency check, not a signature: someone who can
rewrite the whole folder can rewrite the manifest too.

Standard library only. No network.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

MANIFEST = "manifest.json"
VERSION = "version.json"
NCT_ID = re.compile(r"NCT\d{8}")        # use fullmatch: a trailing newline or anything else is not an id
CHUNK = 1 << 16


class SnapshotError(RuntimeError):
    """The snapshot is missing, incomplete or does not match its manifest. Nothing is computed from it."""


_ESCAPES = re.compile(r"\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)|\x1b\[[0-?]*[ -/]*[@-~]|\x1b[P^_X][^\x1b]*\x1b\\|\x1b[@-Z\\-_]")
_CONTROLS = re.compile(r"[\x00-\x1f\x7f-\x9f]")


def clean(value, limit: int = 120) -> str:
    """Text that came from registry pages, a model or a file, made safe to print: terminal escape sequences removed whole, every other
    control character (newline and tab included) turned into a space, and long text cut at `limit` characters."""
    text = _CONTROLS.sub(" ", _ESCAPES.sub("", str(value)))
    return text if len(text) <= limit else text[:limit] + "..."


def valid_nct(value) -> bool:
    return isinstance(value, str) and NCT_ID.fullmatch(value) is not None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_digest(manifest: dict) -> str:
    """sha256 of the manifest without snapshot_sha256, as JSON with sorted keys, no spaces and UTF-8 text."""
    body = {k: v for k, v in manifest.items() if k != "snapshot_sha256"}
    return hashlib.sha256(json.dumps(body, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def _g(d, *path, default=None):
    for p in path:
        if not isinstance(d, dict):
            return default
        d = d.get(p)
        if d is None:
            return default
    return d


def _text(v) -> str:
    return v if isinstance(v, str) else ""


def study_record(study: dict) -> dict:
    """The fields the atlas reads from one API v2 study, under plain names. A missing field becomes None or an empty list, never a guess.
    Primary outcomes keep every entry, in registry order, with measure, description and time frame exactly as written.
    Raises SnapshotError when the NCT id is not of the form NCT followed by 8 digits."""
    ps = (study or {}).get("protocolSection") or {}
    nct = _g(ps, "identificationModule", "nctId")
    if not valid_nct(nct):
        raise SnapshotError(f"a study has no valid NCT id ({str(nct)[:20]!r})")
    outcomes = []
    for i, o in enumerate(_g(ps, "outcomesModule", "primaryOutcomes", default=[]) or [], 1):
        o = o if isinstance(o, dict) else {}
        outcomes.append({"index": i, "measure": _text(o.get("measure")), "description": _text(o.get("description")),
                         "time_frame": _text(o.get("timeFrame"))})
    return {
        "nct_id": nct,
        "title": _g(ps, "identificationModule", "briefTitle"),
        "overall_status": _g(ps, "statusModule", "overallStatus"),
        "start_date": _g(ps, "statusModule", "startDateStruct", "date"),
        "start_type": _g(ps, "statusModule", "startDateStruct", "type"),
        "first_posted": _g(ps, "statusModule", "studyFirstPostDateStruct", "date"),
        "study_type": _g(ps, "designModule", "studyType"),
        "phases": list(_g(ps, "designModule", "phases", default=[]) or []),
        "conditions": [c for c in (_g(ps, "conditionsModule", "conditions", default=[]) or []) if isinstance(c, str)],
        "keywords": [k for k in (_g(ps, "conditionsModule", "keywords", default=[]) or []) if isinstance(k, str)],
        "sponsor": _g(ps, "sponsorCollaboratorsModule", "leadSponsor", "name"),
        "has_results": bool(study.get("hasResults")),
        "primary_outcomes": outcomes,
    }


def outcome_entries(record: dict) -> list[dict]:
    """The tagged unit: one entry per primary outcome, id 'NCT...:P<n>'. A study with no primary outcome gets one empty entry ':P0'
    marked no_primary_outcome, so it is counted as 'not stated' instead of disappearing."""
    nct = record["nct_id"]
    if not record["primary_outcomes"]:
        return [{"entry_id": f"{nct}:P0", "nct_id": nct, "index": 0, "measure": "", "description": "", "time_frame": "",
                 "no_primary_outcome": True}]
    return [{"entry_id": f"{nct}:P{o['index']}", "nct_id": nct, "index": o["index"], "measure": o["measure"],
             "description": o["description"], "time_frame": o["time_frame"], "no_primary_outcome": False}
            for o in record["primary_outcomes"]]


@dataclass
class Snapshot:
    path: Path
    manifest: dict
    version: dict
    pages: dict[str, list[dict]] = field(default_factory=dict)   # route -> pages in order

    @property
    def digest(self) -> str:
        return self.manifest["snapshot_sha256"]

    def routes(self) -> list[str]:
        return list(self.pages)

    def studies(self, route: str) -> list[dict]:
        if route not in self.pages:
            raise SnapshotError(f"the snapshot has no route '{route}' (it has: {', '.join(self.pages) or 'none'})")
        return [s for page in self.pages[route] for s in page.get("studies", [])]

    def records(self, route: str) -> list[dict]:
        return [study_record(s) for s in self.studies(route)]


def _inside(root: Path, rel) -> Path:
    """The file a manifest path names, or SnapshotError when it is not a plain relative path inside the folder."""
    if not isinstance(rel, str) or not rel or "\\" in rel or ":" in rel:
        raise SnapshotError(f"manifest path {str(rel)[:60]!r} is not a plain relative path")
    pp = PurePosixPath(rel)
    if pp.is_absolute() or ".." in pp.parts or rel.startswith("/"):
        raise SnapshotError(f"manifest path {rel[:60]!r} points outside the snapshot")
    target = (root / pp).resolve()
    if not target.is_relative_to(root.resolve()):
        raise SnapshotError(f"manifest path {rel[:60]!r} resolves outside the snapshot")
    return target


def load(path: Path | str, *, verify: bool = True) -> Snapshot:
    """Load a snapshot folder. With verify (the default) every rule in the module docstring is checked first; SnapshotError otherwise."""
    root = Path(path)
    mpath = root / MANIFEST
    if not mpath.is_file():
        raise SnapshotError(f"{root.name}: no {MANIFEST}; a fetch that did not finish leaves none, and nothing is read from it")
    try:
        manifest = json.loads(mpath.read_text(encoding="utf-8"))
        files = manifest["files"]
        routes = manifest["routes"]
        if not isinstance(files, dict) or not isinstance(routes, list):
            raise TypeError("files must be an object and routes a list")
        page_paths = [p["file"] for r in routes for p in r["pages"]]
        route_names = [r["name"] for r in routes]
    except (ValueError, KeyError, TypeError) as exc:
        raise SnapshotError(f"{MANIFEST} could not be read ({type(exc).__name__}: {exc})") from exc
    if verify:
        if manifest_digest(manifest) != manifest.get("snapshot_sha256"):
            raise SnapshotError("the snapshot sha256 does not match the manifest (an edited manifest, or a different snapshot)")
        if len(set(route_names)) != len(route_names):
            raise SnapshotError("a route name is listed twice in the manifest")
        if len(set(page_paths)) != len(page_paths):
            raise SnapshotError("a page is listed twice in the manifest's routes")
        for rel in [VERSION] + page_paths:
            if rel not in files:
                raise SnapshotError(f"{str(rel)[:60]!r} is read but not listed in the manifest's files")
        for rel in files:
            _inside(root, rel)
        # The listed names must be exactly the files on disk, spelled as on disk: this also catches aliases of one file such as
        # 'condition/./page-0001.json', a different letter case or a trailing dot, which the file system may treat as the same file.
        on_disk = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p.name != MANIFEST}
        extra, absent = sorted(on_disk - set(files)), sorted(set(files) - on_disk)
        if extra:
            raise SnapshotError(f"file(s) not named in the manifest: {', '.join(extra[:5])}")
        if absent:
            raise SnapshotError(f"manifest name(s) that are not a file on disk as spelled: {', '.join(str(a)[:60] for a in absent[:5])}")
        for rel, want in files.items():
            fp = root / rel
            if not fp.is_file():
                raise SnapshotError(f"{rel} is named in the manifest but missing")
            if sha256_file(fp) != want:
                raise SnapshotError(f"{rel} does not match its sha256 in the manifest")
    try:
        version = json.loads(_inside(root, VERSION).read_text(encoding="utf-8"))
        pages = {}
        for name, r in zip(route_names, routes):
            pages[name] = [json.loads(_inside(root, p["file"]).read_text(encoding="utf-8")) for p in r["pages"]]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise SnapshotError(f"a snapshot file could not be read ({type(exc).__name__}: {exc})") from exc
    if len(pages) != len(routes):
        raise SnapshotError("a route name is listed twice in the manifest")
    for r in routes:
        name, plist, seen, counted = r["name"], pages[r["name"]], set(), 0
        for meta, page in zip(r["pages"], plist):
            studies = page.get("studies") if isinstance(page, dict) else None
            if not isinstance(studies, list):
                raise SnapshotError(f"route {name}: a page has no studies list")
            if meta.get("studies") != len(studies):
                raise SnapshotError(f"route {name}: {str(meta.get('file'))[:60]} holds {len(studies)} studies, "
                                    f"the manifest says {meta.get('studies')}")
            for s in studies:
                nct = _g(s, "protocolSection", "identificationModule", "nctId")
                if not valid_nct(nct):
                    raise SnapshotError(f"route {name}: a study has no valid NCT id ({str(nct)[:20]!r})")
                if nct in seen:
                    raise SnapshotError(f"route {name}: {nct} appears twice")
                seen.add(nct)
            counted += len(studies)
        if not (counted == r.get("studies_received") == r.get("total_count")):
            raise SnapshotError(f"route {name}: {counted} studies on its pages, but the manifest says {r.get('studies_received')} "
                                f"received and a total of {r.get('total_count')}")
    return Snapshot(path=root, manifest=manifest, version=version, pages=pages)
