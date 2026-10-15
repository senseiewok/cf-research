#!/usr/bin/env python3
"""A reusable intake tool so local agents can propose Resources rows from any permitted site.

Takes a host whose catalog entry says `access: fetch` and a list of addresses.
Conduct rules:
  - Reads robots.txt first (every agent group; unreadable means stop).
  - Loads one page at a time under the project user agent, same-host only, 3 seconds apart, with a size cap.
  - Keeps screenshots and text only in a scratch folder outside every repository and deletes them on request.
  - Asks the local worker (inside a strict data boundary, no tools) for a draft row.
  - Drops any draft that breaks the checker's rules (word cap, no digits or contact details, quotes, sentence breaks).
  - Writes `lead_only` rows. A person confirms; the tool never writes `confirmed_by`.

Dry run on one real host:
  python tools/resources/propose_resources.py --source-id cftrust-resource-pages \\
      --urls https://www.cysticfibrosis.org.uk/what-is-cystic-fibrosis/faqs --dry-run

In dry-run mode, the tool:
  1. Validates the source `cftrust-resource-pages` against `sources/catalog.yaml` (`access: fetch`).
  2. Confirms that all URLs share the landing-page host (`www.cysticfibrosis.org.uk`).
  3. Verifies the scratch directory sits outside all repositories.
  4. Plans the sequence: robots.txt check, paced page fetch, data boundary prompt, and lead_only row generation.
  5. Emits the plan without opening sockets or invoking network calls.
"""

from __future__ import annotations

import argparse
import datetime
import html
import json
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urljoin, urlparse, urlsplit

import yaml

# Resolve roots
REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = REPO_ROOT / "sources" / "catalog.yaml"

# Import checker rules from check_resources if available
try:
    from tools.resources import check_resources
except ImportError:
    # If running from tools/resources directly
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import check_resources  # type: ignore

PROJECT_UA = "senseiewok-research-evidence/1.0.0 (+https://github.com/senseiewok/cf-skills)"
MAX_PAGE_BYTES = 500_000
MAX_ROBOTS_BYTES = 500_000
MAX_REDIRECTS = 3
DEFAULT_PAUSE_S = 3.0
ROBOTS_TIMEOUT_S = 15.0
PAGE_TIMEOUT_S = 30.0

# AI and crawler agent tokens checked in robots.txt
AI_AGENT_TOKENS = frozenset({
    "*",
    "senseiewok-research-evidence",
    "anthropic-ai",
    "claudebot",
    "google-extended",
    "gptbot",
    "chatgpt-user",
    "ccbot",
    "bytespider",
    "cohere-ai",
    "perplexitybot",
    "diffbot",
})


class Refused(RuntimeError):
    """An automated request this tool will not make. Never retried."""


class ConductViolation(RuntimeError):
    """A violation of lab safety, conduct, or repository privacy invariants."""


def is_inside_repo(path: Path) -> bool:
    """True if path resolves inside any repository root in the workspace."""
    resolved = path.resolve()
    # Check parents for .git or known repo names
    curr = resolved
    for _ in range(len(resolved.parts)):
        if (curr / ".git").is_dir() or (curr / ".git").is_file():
            return True
        curr = curr.parent
    return False


def get_default_scratch_dir() -> Path:
    """Returns a scratch directory outside every repository."""
    # Preferred: sibling cf-lab-files/scratch/resource_intake
    env_lab_files = os.environ.get("LAB_FILES")
    if env_lab_files:
        cand = Path(env_lab_files) / "scratch" / "resource_intake"
        if not is_inside_repo(cand):
            return cand

    # Check sibling ../cf-lab-files relative to repository root
    # Note: in a worktree, REPO_ROOT may be .worktrees/cf-research/<branch>
    curr = REPO_ROOT
    for _ in range(5):
        sibling = curr.parent / "cf-lab-files"
        if sibling.is_dir() and not is_inside_repo(sibling):
            return sibling / "scratch" / "resource_intake"
        curr = curr.parent

    # Fallback to system temp directory outside git repos
    import tempfile
    temp_dir = Path(tempfile.gettempdir()) / "cf-lab-scratch" / "resource_intake"
    if is_inside_repo(temp_dir):
        raise ConductViolation(f"Temp directory {temp_dir} resolves inside a repository!")
    return temp_dir


def load_catalog(path: Path | None = None) -> list[dict]:
    p = path or CATALOG_PATH
    if not p.is_file():
        raise Refused(f"catalog not found: {p}")
    with p.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict) or "sources" not in data:
        raise Refused(f"catalog at {p} has invalid structure")
    return [s for s in data["sources"] if isinstance(s, dict)]


def find_catalog_entry(catalog: list[dict], host: str | None = None, source_id: str | None = None) -> dict:
    if source_id:
        for s in catalog:
            if s.get("id") == source_id:
                if s.get("access") != "fetch":
                    raise Refused(f"catalog entry '{source_id}' has access '{s.get('access')}'; fetch refused")
                return s
        raise Refused(f"source '{source_id}' not found in catalog")

    if host:
        host_lower = host.lower()
        matching = []
        for s in catalog:
            entry_host = urlsplit(s.get("landing_page") or s.get("url") or s.get("base_url") or "").hostname
            if entry_host and entry_host.lower() == host_lower:
                matching.append(s)

        if not matching:
            raise Refused(f"no catalog entry found for host '{host}'")

        # Prioritize entries with kind: resource_pages and access: fetch
        resource_fetch = [s for s in matching if s.get("kind") == "resource_pages" and s.get("access") == "fetch"]
        if resource_fetch:
            return resource_fetch[0]

        # Next check any entry with access: fetch
        any_fetch = [s for s in matching if s.get("access") == "fetch"]
        if any_fetch:
            return any_fetch[0]

        # None permitted fetch
        access_list = ", ".join(f"{s.get('id')}:{s.get('access')}" for s in matching)
        raise Refused(f"catalog access for host '{host}' does not permit fetch ({access_list}); automated fetch refused")

    raise Refused("either host or source_id must be provided")


def extract_host(entry: dict) -> str:
    url = entry.get("landing_page") or entry.get("url") or entry.get("base_url") or ""
    host = urlsplit(url).hostname
    if not host:
        raise Refused(f"catalog entry '{entry.get('id')}' has no parseable host")
    return host.lower()


def parse_robots_rules(body: str) -> list[tuple[set[str], list[tuple[str, str]]]]:
    """RFC 9309 parser returning a list of (user_agents, [(allow/disallow, pattern)])."""
    groups: list[tuple[set[str], list[tuple[str, str]]]] = []
    agents: set[str] = set()
    rules: list[tuple[str, str]] = []
    for raw in body.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        key, _, val = line.partition(":")
        key, val = key.strip().lower(), val.strip()
        if key == "user-agent":
            if rules:
                groups.append((agents, rules))
                agents, rules = set(), []
            if val:
                agents.add(val.lower())
        elif key in ("allow", "disallow"):
            rules.append((key, val))
    if agents or rules:
        groups.append((agents, rules))
    return groups


def _path_matches(path: str, pattern: str) -> bool:
    anchored = pattern.endswith("$")
    pat = pattern[:-1] if anchored else pattern
    parts = pat.split("*")
    pos = 0
    for i, part in enumerate(parts):
        if i == 0:
            if not path.startswith(part):
                return False
            pos = len(part)
        else:
            idx = path.find(part, pos)
            if idx < 0:
                return False
            pos = idx + len(part)
    return (pos == len(path)) if anchored else True


def check_agent_rule(rules: list[tuple[str, str]], path: str) -> bool:
    best: tuple[int, str] | None = None
    for key, pattern in rules:
        if pattern == "" or not _path_matches(path, pattern):
            continue
        if best is None or len(pattern) > best[0] or (len(pattern) == best[0] and key == "allow"):
            best = (len(pattern), key)
    return best is None or best[1] == "allow"


def check_all_robots_groups(body: str, path: str) -> tuple[bool, str]:
    """Evaluates robots.txt against project UA, wildcard, and every declared AI agent group."""
    groups = parse_robots_rules(body)
    # Check our own UA token
    our_token = PROJECT_UA.split("/")[0].lower()
    our_rules = [r for a, r in groups if our_token in a]
    if our_rules:
        all_rules = [rule for r in our_rules for rule in r]
        if not check_agent_rule(all_rules, path):
            return False, f"robots.txt disallows {path} for {our_token}"

    # Check wildcard
    star_rules = [r for a, r in groups if "*" in a]
    if star_rules:
        all_rules = [rule for r in star_rules for rule in r]
        if not check_agent_rule(all_rules, path):
            return False, f"robots.txt wildcard disallows {path}"

    # Check every other AI agent group mentioned in the file
    for a, r in groups:
        for token in AI_AGENT_TOKENS:
            if token in a:
                if not check_agent_rule(r, path):
                    return False, f"robots.txt group {a} disallows {path}"
    return True, ""


@dataclass
class FetchResult:
    url: str
    status_code: int
    text: str
    title: str
    bytes_count: int


class IntakeSession:
    """Handles polite HTTP requests with strict conduct rules."""

    def __init__(
        self,
        host: str,
        pause_s: float = DEFAULT_PAUSE_S,
        max_bytes: int = MAX_PAGE_BYTES,
        dry_run: bool = False,
        session=None,
    ):
        self.host = host.lower()
        self.pause_s = pause_s
        self.max_bytes = max_bytes
        self.dry_run = dry_run
        self.session = session
        self._last_call: float | None = None
        self._robots_body: str | None = None

    def _pace(self) -> None:
        if self._last_call is not None:
            elapsed = time.monotonic() - self._last_call
            remaining = self.pause_s - elapsed
            if remaining > 0:
                time.sleep(remaining)
        self._last_call = time.monotonic()

    def load_robots(self) -> str:
        if self._robots_body is not None:
            return self._robots_body
        robots_url = f"https://{self.host}/robots.txt"
        if self.dry_run:
            self._robots_body = ""
            return ""

        import requests

        sess = self.session or requests.Session()
        current = robots_url
        for _ in range(MAX_REDIRECTS + 1):
            self._pace()
            try:
                resp = sess.get(
                    current,
                    headers={"User-Agent": PROJECT_UA},
                    timeout=ROBOTS_TIMEOUT_S,
                    allow_redirects=False,
                )
            except Exception as e:
                raise Refused(f"robots.txt for {self.host} unreadable ({e}); stopping run") from e

            if 300 <= resp.status_code < 400:
                loc = resp.headers.get("Location", "")
                parts = urlsplit(loc)
                if parts.hostname and parts.hostname.lower() != self.host:
                    raise Refused(f"robots.txt for {self.host} redirected cross-host to {loc}; unreadable, stopping")
                if parts.scheme and parts.scheme != "https":
                    raise Refused(f"robots.txt for {self.host} redirected insecurely to {loc}; unreadable, stopping")
                current = urljoin(current, loc)
                continue
            break
        else:
            raise Refused(f"robots.txt for {self.host}: too many redirects; stopping")

        if resp.status_code in (404, 410):
            self._robots_body = ""
            return ""
        if resp.status_code == 200:
            if len(resp.content) > MAX_ROBOTS_BYTES:
                raise Refused(f"robots.txt for {self.host} exceeds {MAX_ROBOTS_BYTES} bytes; unreadable, stopping")
            self._robots_body = resp.text
            return self._robots_body

        # Any other code (401, 403, 500, 502, 503) means unreadable -> STOP!
        raise Refused(f"robots.txt for {self.host} unreadable (status {resp.status_code}); stopping run")

    def fetch_page(self, url: str) -> FetchResult:
        parts = urlsplit(url)
        if parts.scheme != "https":
            raise Refused(f"URL {url} must use https")
        if (parts.hostname or "").lower() != self.host:
            raise Refused(f"URL host '{parts.hostname}' does not match catalog host '{self.host}'")

        # Robots check
        robots = self.load_robots()
        allowed, reason = check_all_robots_groups(robots, parts.path or "/")
        if not allowed:
            raise Refused(f"fetch disallowed: {reason}")

        if self.dry_run:
            return FetchResult(url, 200, f"[Dry Run Content for {url}]", "Dry Run Page", 100)

        import requests

        sess = self.session or requests.Session()
        current = url

        for _ in range(MAX_REDIRECTS + 1):
            self._pace()
            try:
                resp = sess.get(
                    current,
                    headers={"User-Agent": PROJECT_UA},
                    timeout=PAGE_TIMEOUT_S,
                    stream=True,
                    allow_redirects=False,
                )
            except Exception as e:
                raise Refused(f"request to {current} failed ({e})") from e

            if 300 <= resp.status_code < 400:
                loc = resp.headers.get("Location", "")
                target_parts = urlsplit(loc)
                if target_parts.hostname and target_parts.hostname.lower() != self.host:
                    raise Refused(f"cross-host redirect from {current} to {loc} is refused")
                if target_parts.scheme and target_parts.scheme != "https":
                    raise Refused(f"insecure redirect from {current} to {loc} is refused")
                # Re-check robots.txt rules against redirect target path
                target_path = target_parts.path or "/"
                allowed, reason = check_all_robots_groups(robots, target_path)
                if not allowed:
                    raise Refused(f"redirect target {loc} disallowed by robots.txt: {reason}")
                current = urljoin(current, loc)
                continue
            break
        else:
            raise Refused(f"too many redirects fetching {url}")

        if resp.status_code != 200:
            raise Refused(f"server answered {resp.status_code} for {current}; not retried")

        # Size cap check
        cl = resp.headers.get("Content-Length")
        if cl and cl.isdigit() and int(cl) > self.max_bytes:
            raise Refused(f"page {current} declared Content-Length {cl} exceeds cap of {self.max_bytes} bytes")

        chunks = []
        total = 0
        try:
            for chunk in resp.iter_content(chunk_size=16384):
                total += len(chunk)
                if total > self.max_bytes:
                    raise Refused(f"page {current} body exceeded size cap of {self.max_bytes} bytes")
                chunks.append(chunk)
        finally:
            resp.close()

        raw_bytes = b"".join(chunks)
        raw_text = raw_bytes.decode("utf-8", errors="replace")

        # Extract title and clean text
        title_match = re.search(r"<title[^>]*>(.*?)</title>", raw_text, re.IGNORECASE | re.DOTALL)
        title = html.unescape(title_match.group(1).strip()) if title_match else ""

        # Strip scripts, styles, tags
        cleaned = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", raw_text, flags=re.IGNORECASE | re.DOTALL)
        cleaned = re.sub(r"<[^>]+>", " ", cleaned)
        cleaned = html.unescape(cleaned)
        cleaned_text = " ".join(cleaned.split())

        return FetchResult(current, 200, cleaned_text, title, total)


def save_to_scratch(
    scratch_dir: Path,
    slug: str,
    text: str,
    metadata: dict,
) -> tuple[Path, Path]:
    """Persists extracted text and metadata inside scratch folder outside repositories."""
    if is_inside_repo(scratch_dir):
        raise ConductViolation(f"Scratch dir {scratch_dir} cannot be inside a repository!")
    scratch_dir.mkdir(parents=True, exist_ok=True)
    text_file = scratch_dir / f"{slug}.txt"
    meta_file = scratch_dir / f"{slug}.meta.json"
    text_file.write_text(text, encoding="utf-8")
    meta_file.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return text_file, meta_file


def clean_scratch(scratch_dir: Path) -> int:
    """Deletes scratch files on request."""
    if is_inside_repo(scratch_dir):
        raise ConductViolation(f"Refusing to delete {scratch_dir}: resolves inside a repository!")
    if not scratch_dir.exists():
        return 0
    count = 0
    for item in scratch_dir.glob("*"):
        if item.is_file():
            item.unlink()
            count += 1
    if not any(scratch_dir.iterdir()):
        scratch_dir.rmdir()
    return count


def build_worker_prompt(url: str, title: str, text: str) -> str:
    """Builds a data-bounded prompt for the local worker model without external tools."""
    # Truncate text to 1000 words inside boundary
    words = text.split()[:1000]
    excerpt = " ".join(words)

    return f"""You are drafting a candidate resource row for Sensei Ewok's CF Lab Resources register.
You are operating inside a strict data boundary. You have NO tools and NO network access.
Extract information ONLY from the provided text excerpt below. Never extrapolate.

=== BEGIN DATA BOUNDARY ===
Page URL: {url}
Page Title: {title}
Page Text Excerpt:
{excerpt}
=== END DATA BOUNDARY ===

Produce a valid JSON object with the following fields:
- "name": string, the concise official name of the resource/service as written by the organisation.
- "summary": string, exactly ONE sentence in the lab's own words describing the support provided. Rules:
    * At most 25 words.
    * Exactly one sentence (no newline, no period except at end).
    * NO double quotes or curly quotes.
    * NO digits except a four-digit year (e.g. 2026).
    * NO URL-like text (http, www, .org, .com).
    * NO email addresses or phone numbers.
- "audience": non-empty list from: [adult_with_cf, teen_young_adult, parent_carer_of_child, partner_family_carer_of_adult, care_team, researcher, everyone]
- "category": non-empty list from: [urgent_crisis, newly_diagnosed_basics, find_care_team, mental_health, medicines_access, money_insurance, daily_life, growing_up_transition, relationships_family_planning, transplant_advanced, grief_bereavement, community_peer, trials_research_registries, advocacy_involved]
- "who_for": exactly one of [stated_anyone, stated_countries, stated_region, not_stated, not_checked]
- "format": exactly one of [online, in_person, both]
- "countries": non-empty list of ISO 3166-1 alpha-2 codes (e.g. ["GB"]) or region values ["europe", "global"]
- "lang": BCP-47 language code (e.g. "en")
- "organisation": name of the publisher organisation

Output JSON only. No explanations.
"""


def validate_draft(draft: dict, source_id: str, catalog_entry: dict) -> tuple[bool, str]:
    """Validates draft row against check_resources rules. Returns (is_valid, reason)."""
    if not isinstance(draft, dict):
        return False, "draft is not a dictionary"

    # Name check
    name = draft.get("name")
    if not isinstance(name, str) or not name.strip():
        return False, "name is missing or empty"
    if check_resources.has_contact_detail(name):
        return False, "name contains phone or email details"

    # Summary check
    summary = draft.get("summary")
    summary_errors = check_resources.check_summary(summary, "summary")
    if summary_errors:
        return False, f"summary invalid: {'; '.join(summary_errors)}"
    if check_resources.has_contact_detail(summary):
        return False, "summary contains phone or email details"

    # Enum checks
    audience = draft.get("audience")
    if not isinstance(audience, list) or not audience:
        return False, "audience must be a non-empty list"
    for a in audience:
        if a not in check_resources.AUDIENCES:
            return False, f"audience '{a}' is not in allowed list"

    category = draft.get("category")
    if not isinstance(category, list) or not category:
        return False, "category must be a non-empty list"
    for c in category:
        if c not in check_resources.CATEGORIES:
            return False, f"category '{c}' is not in allowed list"

    who_for = draft.get("who_for")
    if who_for not in check_resources.WHO_FOR:
        return False, f"who_for '{who_for}' is not allowed"

    fmt = draft.get("format", "online")
    if fmt not in check_resources.FORMATS:
        return False, f"format '{fmt}' is not allowed"

    countries = draft.get("countries")
    if not isinstance(countries, list) or not countries:
        return False, "countries must be a non-empty list"
    for code in countries:
        if code not in check_resources.ISO_COUNTRIES and code not in check_resources.REGIONS:
            return False, f"country '{code}' is not an ISO code or region"

    return True, ""


def slugify(text: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9]+", "-", text).strip("-").lower()
    return cleaned or "resource"


def build_lead_only_row(
    draft: dict,
    url: str,
    source_id: str,
    today: datetime.date | None = None,
) -> dict:
    """Constructs a validated lead_only row. Strictly forbids confirmed_by."""
    check_date = today or datetime.date.today()
    rid = slugify(draft["name"])

    # INVARIANT: A person confirms; the tool never writes confirmed_by or content_hash.
    forbidden = ("confirmed_by", "confirmed_on", "content_hash", "summary_attested_by", "summary_attested_on")
    for k in forbidden:
        if k in draft:
            raise ConductViolation(f"The intake tool writes lead_only rows and strictly never writes {k}!")

    row = {
        "id": rid,
        "name": draft["name"].strip(),
        "organisation": draft.get("organisation", "").strip(),
        "countries": draft["countries"],
        "url": url,
        "audience": draft["audience"],
        "category": draft["category"],
        "summary": draft["summary"].strip(),
        "who_for": draft["who_for"],
        "format": draft.get("format", "online"),
        "lang": draft.get("lang", "en"),
        "source_id": source_id,
        "basis": "lead_only",
        "last_checked": check_date.isoformat(),
        "status": "active",
    }
    return row


def propose_resources(
    urls: list[str],
    source_id: str | None = None,
    host: str | None = None,
    catalog_path: Path | None = None,
    scratch_dir: Path | None = None,
    pause_s: float = DEFAULT_PAUSE_S,
    max_bytes: int = MAX_PAGE_BYTES,
    dry_run: bool = False,
    worker_fn=None,
    delete_scratch_after: bool = False,
    session=None,
) -> list[dict]:
    """Orchestrates intake of resources from permitted sites into lead_only rows."""
    catalog = load_catalog(catalog_path)
    entry = find_catalog_entry(catalog, host=host, source_id=source_id)
    canonical_host = extract_host(entry)
    sid = entry["id"]

    scratch = scratch_dir or get_default_scratch_dir()
    if is_inside_repo(scratch):
        raise ConductViolation(f"Scratch directory {scratch} is inside a git repo!")

    intake_session = IntakeSession(
        host=canonical_host,
        pause_s=pause_s,
        max_bytes=max_bytes,
        dry_run=dry_run,
        session=session,
    )

    valid_rows = []

    try:
        for url in urls:
            # Load and parse page
            fetch_res = intake_session.fetch_page(url)
            slug = slugify(fetch_res.title or urlsplit(url).path)

            # Save text and metadata to scratch outside repository
            if not dry_run:
                meta = {
                    "url": url,
                    "title": fetch_res.title,
                    "bytes": fetch_res.bytes_count,
                    "host": canonical_host,
                    "source_id": sid,
                    "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                }
                save_to_scratch(scratch, slug, fetch_res.text, meta)

            # Prompt worker
            prompt = build_worker_prompt(url, fetch_res.title, fetch_res.text)

            if dry_run:
                draft = {
                    "name": f"Help Service for {fetch_res.title or 'CF'}",
                    "summary": "Information and support resources provided for individuals and families affected by cystic fibrosis.",
                    "audience": ["adult_with_cf", "parent_carer_of_child"],
                    "category": ["daily_life"],
                    "who_for": "stated_anyone",
                    "format": "online",
                    "countries": ["GB"],
                    "lang": "en",
                    "organisation": entry.get("publisher", "CF Organisation"),
                }
            elif worker_fn is not None:
                raw_response = worker_fn(prompt)
                if isinstance(raw_response, str):
                    # Clean possible markdown block
                    cleaned = re.sub(r"^```(?:json)?\s*", "", raw_response.strip())
                    cleaned = re.sub(r"\s*```$", "", cleaned)
                    try:
                        draft = json.loads(cleaned)
                    except json.JSONDecodeError:
                        print(f"WARN: Dropping draft for {url}: Worker output was not valid JSON", file=sys.stderr)
                        continue
                else:
                    draft = raw_response
            else:
                raise Refused("Worker callable not configured and live model not specified")

            # Validate draft against checker rules
            is_valid, reason = validate_draft(draft, sid, entry)
            if not is_valid:
                print(f"WARN: Dropping draft for {url}: {reason}", file=sys.stderr)
                continue

            # Form lead_only row
            row = build_lead_only_row(draft, url, sid)
            valid_rows.append(row)
    finally:
        if delete_scratch_after and not dry_run:
            clean_scratch(scratch)

    return valid_rows


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="A reusable intake tool so local agents can propose Resources rows from permitted sites."
    )
    parser.add_argument("--urls", nargs="*", default=[], help="List of addresses to propose")
    parser.add_argument("--source-id", help="Source ID in catalog.yaml (must have access: fetch)")
    parser.add_argument("--host", help="Host name matching a catalog entry with access: fetch")
    parser.add_argument("--catalog", type=Path, default=CATALOG_PATH, help="Path to catalog.yaml")
    parser.add_argument("--scratch-dir", type=Path, help="Override scratch directory (must be outside all repositories)")
    parser.add_argument("--pause", type=float, default=DEFAULT_PAUSE_S, help="Pause in seconds between requests (default 3.0)")
    parser.add_argument("--max-bytes", type=int, default=MAX_PAGE_BYTES, help="Size cap per page (default 500,000)")
    parser.add_argument("--output", type=Path, help="Write proposed YAML rows to this file")
    parser.add_argument("--dry-run", action="store_true", help="Simulate execution without opening sockets")
    parser.add_argument("--delete-scratch", action="store_true", help="Delete scratch files after execution")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.delete_scratch and not args.urls:
        scratch = args.scratch_dir or get_default_scratch_dir()
        count = clean_scratch(scratch)
        print(f"Cleaned {count} files from scratch directory {scratch}")
        return 0

    if not args.urls:
        print("ERROR: the following arguments are required: --urls (unless using --delete-scratch)", file=sys.stderr)
        return 1
    try:
        rows = propose_resources(
            urls=args.urls,
            source_id=args.source_id,
            host=args.host,
            catalog_path=args.catalog,
            scratch_dir=args.scratch_dir,
            pause_s=args.pause,
            max_bytes=args.max_bytes,
            dry_run=args.dry_run,
            delete_scratch_after=args.delete_scratch,
        )
    except (Refused, ConductViolation) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    yaml_text = yaml.dump(rows, sort_keys=False, allow_unicode=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as fh:
            fh.write(yaml_text)
        print(f"Wrote {len(rows)} lead_only rows to {args.output}")
    else:
        print(yaml_text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
