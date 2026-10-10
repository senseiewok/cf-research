"""Read a trial-atlas snapshot written by fetch_snapshot.py, check it against its manifest, and turn its studies into plain records.

A snapshot is a folder:  manifest.json, version.json and one folder per retrieval route (condition/, term/) of page-NNNN.json files.
The manifest holds the sha256 of every other file and the sha256 of the whole snapshot (see snapshot_digest). Loading checks every hash
first and refuses a snapshot that does not match: a count is only ever computed from the bytes the manifest names.

Standard library only. No network.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

MANIFEST = "manifest.json"
VERSION = "version.json"
NCT_ID = re.compile(r"^NCT\d{8}$")
CHUNK = 1 << 16


class SnapshotError(RuntimeError):
    """The snapshot is missing, incomplete or does not match its manifest. Nothing is computed from it."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def snapshot_digest(file_hashes: dict[str, str]) -> str:
    """One hash for the whole snapshot: sha256 of the lines 'relative/path<TAB>sha256' sorted by path, each ending in a newline."""
    lines = "".join(f"{p}\t{file_hashes[p]}\n" for p in sorted(file_hashes))
    return hashlib.sha256(lines.encode("utf-8")).hexdigest()


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
    Primary outcomes keep every entry, in registry order, with measure, description and time frame exactly as written."""
    ps = study.get("protocolSection") or {}
    outcomes = []
    for i, o in enumerate(_g(ps, "outcomesModule", "primaryOutcomes", default=[]) or [], 1):
        o = o if isinstance(o, dict) else {}
        outcomes.append({"index": i, "measure": _text(o.get("measure")), "description": _text(o.get("description")),
                         "time_frame": _text(o.get("timeFrame"))})
    return {
        "nct_id": _g(ps, "identificationModule", "nctId"),
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


def load(path: Path | str, *, verify: bool = True) -> Snapshot:
    """Load a snapshot folder. With verify (the default) every file is re-hashed and must match the manifest, the snapshot hash must
    match, and no page file may exist that the manifest does not name. Raises SnapshotError otherwise."""
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
    except (ValueError, KeyError, TypeError) as exc:
        raise SnapshotError(f"{MANIFEST} could not be read ({type(exc).__name__}: {exc})") from exc
    if verify:
        on_disk = {p.relative_to(root).as_posix() for p in root.rglob("*.json") if p.name != MANIFEST}
        extra = sorted(on_disk - set(files))
        if extra:
            raise SnapshotError(f"file(s) not named in the manifest: {', '.join(extra[:5])}")
        for rel, want in files.items():
            fp = root / rel
            if not fp.is_file():
                raise SnapshotError(f"{rel} is named in the manifest but missing")
            if sha256_file(fp) != want:
                raise SnapshotError(f"{rel} does not match its sha256 in the manifest")
        if snapshot_digest(files) != manifest.get("snapshot_sha256"):
            raise SnapshotError("the snapshot sha256 does not match the file hashes in the manifest")
    try:
        version = json.loads((root / VERSION).read_text(encoding="utf-8"))
        pages = {}
        for r in routes:
            pages[r["name"]] = [json.loads((root / p["file"]).read_text(encoding="utf-8")) for p in r["pages"]]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise SnapshotError(f"a snapshot file could not be read ({type(exc).__name__}: {exc})") from exc
    return Snapshot(path=root, manifest=manifest, version=version, pages=pages)
