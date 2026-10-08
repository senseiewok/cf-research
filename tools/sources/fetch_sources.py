#!/usr/bin/env python3
"""Fetch the source documents this repo is permitted to fetch, and explain the rest.

Reads sources/catalog.yaml. Downloads entries whose access policy is 'fetch'.
Prints browser instructions for 'manual' entries. Refuses 'forbidden' ones.
Writes sources/manifest.json recording the SHA-256 of every file on disk (and keeping
the hash an earlier run recorded for a file not on this machine), so an extraction
can be reproduced without redistributing the documents themselves.

Usage:
    python fetch_sources.py                 fetch permitted sources, list the rest
    python fetch_sources.py --manual        print browser instructions only
    python fetch_sources.py --verify        re-hash local files against the manifest (add --only ID [ID ...] to verify a subset)
    python fetch_sources.py --record        checksum files already on disk (e.g. saved by hand); no network
    python fetch_sources.py --only ID [ID]  restrict to specific catalog ids
    python fetch_sources.py --dest PATH     write somewhere other than sources/downloads
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin, urlparse

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG = REPO_ROOT / "sources" / "catalog.yaml"
MANIFEST = REPO_ROOT / "sources" / "manifest.json"
DEFAULT_DEST = REPO_ROOT / "sources" / "downloads"

# Only these two policies permit an automated request. See sources/README.md.
AUTOMATABLE = {"fetch", "api"}
CHUNK = 1 << 16
# Download conduct (see Downloader). The caps are generous for a PDF report and tight for anything else.
MAX_DOWNLOAD_BYTES = 200_000_000
MAX_ROBOTS_BYTES = 500_000
MAX_REDIRECTS = 3
PAUSE_S = 2.0
MAX_CRAWL_DELAY_S = 30.0
TIMEOUT_S = 120


def load_catalog(path: Path | None = None) -> dict:
    with (path or CATALOG).open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


class Refused(RuntimeError):
    """A download this tool will not make, with the reason. Never retried."""


def _evidence_dirs() -> list[Path]:
    env = os.environ.get("EVIDENCE_SKILL_DIR")
    return ([Path(env)] if env else []) + [REPO_ROOT.parent / "cf-skills" / ".claude" / "skills" / "cf-evidence-loop" / "scripts"]


def load_evidence():
    """The evidence skill's HTTP module: it holds the one maintained robots rule (AI-agent groups included) and the project user agent."""
    for d in _evidence_dirs():
        if (d / "evidence" / "http.py").is_file():
            sys.path.insert(0, str(d))
            try:
                from evidence import http as module  # noqa: PLC0415
            except ImportError as exc:
                raise Refused(f"the evidence skill could not be imported ({exc}); install this folder's requirements.txt") from exc
            finally:
                sys.path.remove(str(d))
            return module
    raise Refused("the evidence skill that holds the robots rule is not available; clone cf-skills next to this repository or set "
                  "EVIDENCE_SKILL_DIR, or download the file by hand (see --manual)")


def evidence_available() -> bool:
    try:
        load_evidence()
    except Refused:
        return False
    return True


def crawl_delay(body: str) -> float:
    """The largest Crawl-delay in a robots.txt, in seconds; 0 when there is none or it is not a positive number."""
    best = 0.0
    for raw in body.splitlines():
        key, _, val = raw.split("#", 1)[0].partition(":")
        if key.strip().lower() == "crawl-delay":
            try:
                best = max(best, float(val.strip()))
            except ValueError:
                pass
    return best


class Downloader:
    """One polite download at a time, under the evidence skill's network rules: the project user agent and nothing else; the host's
    robots.txt read once and obeyed (a group that names an AI agent applies, an unreadable file refuses); https only; redirects only
    within the same host and checked against robots again; one request per file with no retry; a pause between requests to a host
    (at least PAUSE_S, or its Crawl-delay up to MAX_CRAWL_DELAY_S); a size cap; no partial file left behind."""

    def __init__(self, session=None, sleep=time.sleep, evidence=None, load_evidence=load_evidence, max_bytes: int = MAX_DOWNLOAD_BYTES):
        self.session, self.sleep, self.evidence, self._load_evidence, self.max_bytes = session, sleep, evidence, load_evidence, max_bytes
        self._robots: dict[str, tuple[str, float]] = {}
        self._requests: dict[str, int] = {}

    def _ev(self):
        if self.evidence is None:
            self.evidence = self._load_evidence()
        return self.evidence

    def _get(self, url: str, *, stream: bool, delay: float):
        if self.session is None:
            import requests  # imported here so --manual, --verify and --record work without it
            self.session = requests.Session()
        host = urlparse(url).netloc.lower()
        if self._requests.get(host):
            self.sleep(min(max(PAUSE_S, delay), MAX_CRAWL_DELAY_S))
        self._requests[host] = self._requests.get(host, 0) + 1
        return self.session.get(url, headers={"User-Agent": self._ev().PROJECT_UA}, timeout=TIMEOUT_S, stream=stream, allow_redirects=False)

    @staticmethod
    def _same_host_https(base: str, location: str, what: str) -> str:
        target = urljoin(base, location or "")
        a, b = urlparse(base), urlparse(target)
        if b.scheme != "https" or b.netloc.lower() != a.netloc.lower():
            raise Refused(f"{what} redirect to {target[:80]} refused: only https redirects within the same host are followed")
        return target

    def _load_robots(self, host: str) -> tuple[str, float]:
        if host in self._robots:
            return self._robots[host]
        url = f"https://{host}/robots.txt"
        for _ in range(MAX_REDIRECTS + 1):
            try:
                resp = self._get(url, stream=False, delay=0.0)
            except OSError as exc:
                raise Refused(f"robots.txt for {host} could not be read ({type(exc).__name__}); nothing is fetched from this host") from exc
            if 300 <= resp.status_code < 400:
                url = self._same_host_https(url, resp.headers.get("Location", ""), "robots.txt")
                continue
            break
        else:
            raise Refused(f"robots.txt for {host}: too many redirects")
        if resp.status_code in (404, 410):
            body = ""
        elif resp.status_code == 200:
            if len(resp.content) > MAX_ROBOTS_BYTES:
                raise Refused(f"robots.txt for {host} is larger than {MAX_ROBOTS_BYTES:,} bytes; nothing is fetched from this host")
            body = resp.text
        else:
            raise Refused(f"robots.txt for {host} answered {resp.status_code}; nothing is fetched from this host")
        self._robots[host] = (body, crawl_delay(body))
        return self._robots[host]

    def robots_verdict(self, url: str) -> tuple[bool, str]:
        """Would this tool be allowed to fetch the URL? Reads the host's robots.txt (once) and nothing else."""
        parts = urlparse(url)
        try:
            body, _ = self._load_robots(parts.netloc.lower())
            ev = self._ev()
        except Refused as exc:
            return False, str(exc)
        if ev.robots_allows(body, parts.path or "/", ev.PROJECT_UA.split("/")[0]):
            return True, ""
        return False, f"robots.txt for {parts.netloc} disallows {parts.path} for this tool (groups that name AI agents apply to it)"

    def download(self, url: str, dest: Path) -> None:
        parts = urlparse(url)
        if parts.scheme != "https" or not parts.netloc:
            raise Refused("only https URLs are downloaded")
        ev = self._ev()
        host = parts.netloc.lower()
        body, delay = self._load_robots(host)
        token = ev.PROJECT_UA.split("/")[0]
        current = url
        for _ in range(MAX_REDIRECTS + 1):
            if not ev.robots_allows(body, urlparse(current).path or "/", token):
                raise Refused(f"robots.txt for {host} disallows {urlparse(current).path} for this tool (groups that name AI agents apply to it)")
            resp = self._get(current, stream=True, delay=delay)
            if 300 <= resp.status_code < 400:
                resp.close()
                current = self._same_host_https(current, resp.headers.get("Location", ""), "document")
                continue
            break
        else:
            raise Refused("document: too many redirects")
        try:
            if resp.status_code != 200:
                raise Refused(f"the server answered {resp.status_code}; not retried")
            declared = resp.headers.get("Content-Length", "")
            if declared.isdigit() and int(declared) > self.max_bytes:
                raise Refused(f"the file is larger than {self.max_bytes:,} bytes")
            tmp = dest.with_suffix(dest.suffix + ".part")
            total = 0
            try:
                with tmp.open("wb") as fh:
                    for chunk in resp.iter_content(CHUNK):
                        total += len(chunk)
                        if total > self.max_bytes:
                            raise Refused(f"the file is larger than {self.max_bytes:,} bytes")
                        fh.write(chunk)
                tmp.replace(dest)
            finally:
                tmp.unlink(missing_ok=True)
        finally:
            resp.close()


def audit_robots(entries: list[dict], downloader: Downloader) -> list[tuple[str, str, str]]:
    """For every `fetch` entry with a URL: (id, "allowed" or "REFUSED", the reason). One robots.txt request per host; no document is requested."""
    rows = []
    for e in entries:
        if e.get("access") != "fetch" or not e.get("url"):
            continue
        ok, why = downloader.robots_verdict(e["url"])
        rows.append((e["id"], "allowed" if ok else "REFUSED", why))
    return rows


_DOWNLOADER: Downloader | None = None


def download(url: str, dest: Path) -> None:
    """Download one permitted file with the Downloader's conduct; the Downloader (its robots cache and pacing) lasts for the process."""
    global _DOWNLOADER
    if _DOWNLOADER is None:
        _DOWNLOADER = Downloader()
    _DOWNLOADER.download(url, dest)


def fetchable(entry: dict) -> bool:
    return entry.get("access") in AUTOMATABLE and bool(entry.get("url"))


def manual_instructions(entries: list[dict]) -> str:
    by_publisher: dict[str, list[dict]] = {}
    for e in entries:
        by_publisher.setdefault(e.get("publisher", "unknown"), []).append(e)

    out = ["", "Manual downloads needed", "=" * 23, ""]
    for publisher, items in sorted(by_publisher.items()):
        out.append(f"{publisher}")
        landing = {i.get("landing_page") for i in items if i.get("landing_page")}
        for page in sorted(landing):
            out.append(f"  Start here: {page}")
        out.append("")
        for i in sorted(items, key=lambda x: str(x.get("data_year", ""))):
            flag = "  (!) " if i.get("warning") else "      "
            target = i.get("filename") or "(no file; link only)"
            out.append(f"{flag}{i['id']}  ->  {target}")
            out.append(f"        {i['title']}")
            if i.get("url"):
                out.append(f"        {i['url']}")
            if i.get("warning"):
                out.append(f"        warning: {i['warning']}")
        out.append("")
    out.append("Save each file into the downloads directory under the filename shown, run --record to checksum what you saved, then --verify.")
    return "\n".join(out)


def previous_hashes() -> dict[str, dict]:
    """The hashes the current manifest holds, by id: {id: {filename, sha256, size_bytes}}. Empty when there is no manifest yet.
    Raises ValueError when the manifest exists but cannot be read, so it is never overwritten (and its hashes lost) by mistake."""
    if not MANIFEST.exists():
        return {}
    try:
        files = json.loads(MANIFEST.read_text(encoding="utf-8"))["files"]
        return {r["id"]: r for r in files if r.get("sha256")}
    except (ValueError, KeyError, TypeError) as exc:
        raise ValueError(f"{MANIFEST.name} exists but could not be read ({type(exc).__name__}); fix or remove it, then run again") from exc


def build_manifest(catalog: dict, dest: Path, previous: dict[str, dict] | None = None) -> dict:
    """Hash every catalog file that is on disk. A file that is not on disk keeps the hash an earlier run recorded for the same id and
    filename (marked present: false), so a run on a machine that holds only some of the files does not erase the others' hashes.
    A hash is never made up: an entry that was never recorded has none."""
    previous = previous or {}
    records = []
    for entry in catalog.get("sources", []):
        name = entry.get("filename")
        if not name:
            continue
        path = dest / name
        rec = {
            "id": entry["id"],
            "filename": name,
            "access": entry.get("access"),
            "url": entry.get("url"),
            "present": path.exists(),
        }
        old = previous.get(entry["id"])
        if path.exists():
            rec["sha256"] = sha256(path)
            rec["size_bytes"] = path.stat().st_size
        elif old and old.get("filename") == name:
            rec["sha256"] = old["sha256"]
            if "size_bytes" in old:
                rec["size_bytes"] = old["size_bytes"]
        records.append(rec)
    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "catalog_updated": str(catalog.get("updated")),
        # Never write a machine path into the tracked manifest.
        "downloads_dir": dest.relative_to(REPO_ROOT).as_posix() if dest.is_relative_to(REPO_ROOT) else "(outside the repository)",
        "files": records,
    }


def record(catalog: dict, dest: Path) -> int:
    """Hash whatever is on disk into the manifest, keeping earlier hashes of files not on disk. Never touches the network."""
    try:
        manifest = build_manifest(catalog, dest, previous_hashes())
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    have =[f for f in manifest["files"] if f["present"]]
    kept = [f for f in manifest["files"] if not f["present"] and f.get("sha256")]
    print(f"recorded {len(have)} of {len(manifest['files'])} files in {MANIFEST.name}"
          + (f"; kept the earlier hash of {len(kept)} file(s) not on disk" if kept else ""))
    manual_missing = [f for f in manifest["files"] if f["access"] == "manual" and not f["present"]]
    for f in manual_missing:
        print(f"still needed (save by hand): {f['id']}  ->  {f['filename']}")
    return 0


def verify(dest: Path, only_ids: list[str] | None = None) -> int:
    if not MANIFEST.exists():
        print("No manifest yet. Run without --verify first.")
        return 1
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))

    wanted = set(only_ids) if only_ids is not None else None
    if wanted is not None:
        known_ids = {rec["id"] for rec in manifest["files"]}
        unknown = wanted - known_ids
        if unknown:
            for uid in sorted(unknown):
                print(f"unknown id: {uid}")
            return 2

    bad = missing = ok = nohash = 0
    for rec in manifest["files"]:
        if wanted is not None and rec["id"] not in wanted:
            continue
        path = dest / rec["filename"]
        if not rec.get("sha256"):
            print(f"NO HASH  {rec['id']}  (no recorded hash yet; save the file, then run --record)")
            nohash += 1
            continue
        if not path.exists():
            print(f"MISSING  {rec['id']}  {rec['filename']}")
            missing += 1
        elif sha256(path) != rec["sha256"]:
            print(f"CHANGED  {rec['id']}  {rec['filename']}")
            bad += 1
        else:
            ok += 1

    print(f"\n{ok} match, {bad} changed, {missing} missing" + (f", {nohash} with no hash" if nohash else ""))
    if bad or missing:
        return 1
    if ok == 0:
        print("nothing was verified: no file had a recorded hash to check. Save the files, run --record, then --verify.")
        return 1
    # A named entry with no hash fails (--only); in a full run, entries not yet recorded are listed but do not fail the others.
    if nohash and wanted is not None:
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    ap.add_argument("--only", nargs="+", metavar="ID")
    ap.add_argument("--manual", action="store_true", help="print browser instructions and exit")
    ap.add_argument("--verify", action="store_true", help="re-hash local files against the manifest")
    ap.add_argument("--force", action="store_true", help="re-download files that already exist")
    ap.add_argument("--audit-robots", action="store_true", help="read each fetch host's robots.txt once (no document is requested) and report whether this tool may fetch every fetch entry under the current rule; exit 1 if any is refused")
    ap.add_argument("--record", action="store_true", help="checksum the files already in the downloads directory (for example ones saved by hand) into the manifest; no network")
    args = ap.parse_args()

    args.dest.mkdir(parents=True, exist_ok=True)
    if args.verify:
        return verify(args.dest, args.only)

    catalog = load_catalog()
    if args.record:
        return record(catalog, args.dest)
    entries = catalog.get("sources", [])
    if args.only:
        wanted = set(args.only)
        entries = [e for e in entries if e["id"] in wanted]
        unknown = wanted - {e["id"] for e in entries}
        if unknown:
            print(f"Unknown catalog id(s): {', '.join(sorted(unknown))}", file=sys.stderr)
            return 2

    if args.audit_robots:
        rows = audit_robots(entries, Downloader())
        for sid, verdict, why in rows:
            print(f"{verdict:8} {sid}{'  ' + why if why else ''}")
        refused = sum(1 for r in rows if r[1] == "REFUSED")
        print(f"\n{len(rows) - refused} allowed, {refused} refused")
        return 1 if refused else 0

    manual = [e for e in entries if e.get("access") == "manual"]
    blocked = [e for e in entries if e.get("access") in {"forbidden", "request"}]

    if args.manual:
        print(manual_instructions(manual))
        return 0

    to_fetch = [e for e in entries if fetchable(e) and e.get("filename")]
    failed: list[str] = []
    for entry in to_fetch:
        path = args.dest / entry["filename"]
        if path.exists() and not args.force:
            print(f"have     {entry['id']}")
            continue
        print(f"fetching {entry['id']} ... ", end="", flush=True)
        try:
            download(entry["url"], path)
            print(f"{path.stat().st_size:,} bytes")
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"FAILED ({type(exc).__name__}{': ' + str(exc) if isinstance(exc, Refused) else ''})")
            failed.append(entry["id"])

    try:
        manifest = build_manifest(catalog, args.dest, previous_hashes())
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    present = sum(1 for f in manifest["files"] if f["present"])
    print(f"\nmanifest: {MANIFEST.relative_to(REPO_ROOT)}  ({present} of {len(manifest['files'])} files present)")

    for entry in blocked:
        reason = "publisher refuses automated access" if entry["access"] == "forbidden" else "needs an application or written permission"
        print(f"skipped  {entry['id']}: {reason}")

    if manual:
        print(manual_instructions(manual))
    if failed:
        print(f"\n{len(failed)} download(s) failed: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
