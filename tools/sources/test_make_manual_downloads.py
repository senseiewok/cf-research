"""Test for make_manual_downloads.py (board T-0063): renders a fixture catalog and checks the page, then the command line, then the
real catalog and the committed page. Usage: python test_make_manual_downloads.py [SCRIPT]
One FAIL line per failing check; prints VERIFIED and exits 0 when all pass."""
import importlib.util
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import unittest
import yaml

DEFAULT_SCRIPT = Path(__file__).resolve().parent / "make_manual_downloads.py"

PROMPT = """I have saved the manual downloads into sources/downloads/ under the file names in manual-downloads.md.
1. Run `python tools/sources/fetch_sources.py --record`, then `python tools/sources/fetch_sources.py --verify`, and tell me which files are missing or changed. Do not fetch anything over the network.
2. For each newly recorded file, open its title page and read it. Fill data_year, pages and text_layer in sources/catalog.yaml from what the document says, quoting the title-page text. Set claim_label to verified only for what you read.
3. Run `python tools/sources/check_catalog_fields.py`, then sync the skill's catalog copy with `python tools/sources/sync_skill_catalog.py --source sources/catalog.yaml --target ../cf-skills/.claude/skills/cf-evidence-loop/catalog.yaml`.
4. Show me the diff and wait for my approval before you commit anything."""


def run_checks(script=None):
    """Run every check against SCRIPT (default: the real script next to this file). Returns the FAIL lines; empty when all pass."""
    cand = Path(script).resolve() if script else DEFAULT_SCRIPT
    try:
        spec = importlib.util.spec_from_file_location("cand_md", cand)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    except Exception as e:  # noqa: BLE001
        return [f"FAIL import: {type(e).__name__}: {str(e)[:160]}"]


    LONG = "Very long sentence without a stop " * 12 + "end. Tail sentence that must not appear."
    CAT = {"schema_version": 1, "sources": [
        {"id": "z-one", "access": "manual", "publisher": "Zeta Org", "title": "Report One", "landing_page": "https://zeta.example/reports", "url": "https://zeta.example/files/one.pdf",
         "filename": "ONE.pdf", "notes": "First sentence here. Second sentence stays out."},
        {"id": "a-two", "access": "manual", "publisher": "alpha Foundation", "title": "Pipe | title — café", "landing_page": "https://alpha.example/", "filename": "TWO.pdf",
         "warning": "filename_year_mismatch", "notes": "Line one\nline two. More."},
        {"id": "z-three", "access": "manual", "publisher": "Zeta Org", "title": "Link only", "landing_page": "https://zeta.example/linkonly"},
        {"id": "n-four", "access": "manual", "title": "No publisher", "landing_page": "https://nopub.example/", "filename": "FOUR.pdf", "notes": LONG},
        {"id": "r-one", "access": "request", "publisher": "Registry Hospital", "title": "Ask Registry", "landing_page": "https://registry.example/"},
        {"id": "f-one", "access": "forbidden", "publisher": "Closed DB", "title": "Closed Database", "landing_page": "https://closed.example/"},
        {"id": "fetch-x", "access": "fetch", "publisher": "Open Pub", "title": "Fetchable report", "landing_page": "https://open.example/", "filename": "X.pdf"},
        {"id": "api-y", "access": "api", "publisher": "Api Pub", "title": "An API", "landing_page": "https://api.example/"},
    ]}

    fails = []


    def check(name, ok, detail=""):
        if not ok:
            fails.append(f"FAIL {name} {detail}".strip())


    def cells(row):
        parts = re.split(r"(?<!\\)\|", row.strip())
        return [c.strip() for c in parts[1:-1]]


    try:
        md = mod.render(CAT)
    except Exception as e:  # noqa: BLE001
        return [f"FAIL render raised {type(e).__name__}: {str(e)[:160]}"]
    check("render returns str", isinstance(md, str))
    lines = md.split("\n")
    check("first line", lines[0] == "# Manual downloads", repr(lines[0]))
    check("generated line", "Generated from sources/catalog.yaml by tools/sources/make_manual_downloads.py. Do not edit by hand." in lines)
    check("count line", "4 to download by hand, 1 to ask for, 1 not to download." in lines)
    check("ends with exactly one newline", md.endswith("\n") and not md.endswith("\n\n"))
    check("no carriage returns", "\r" not in md)
    check("no trailing spaces", not any(l.endswith(" ") for l in lines))
    check("deterministic", mod.render(CAT) == md)
    for h in ("## How this works", "## Download in a browser", "## Ask first (access by request)", "## Do not download", "## When you have saved the files"):
        check(f"heading {h}", h in lines)
    check("how-this-works steps", "1. Open the link for a source in your browser." in lines and "2. Save the file into `sources/downloads/` under the exact file name shown." in lines
          and "3. When you have saved the files, paste the prompt at the bottom into your agent." in lines)

    hs = [l for l in lines if l.startswith("### ")]
    check("publisher headings sorted case-insensitively, unknown publisher named", hs == ["### alpha Foundation", "### Unknown publisher", "### Zeta Org"], str(hs))
    check("table header", lines.count("| ID | What | Open | Save as | Notes |") == 3 and lines.count("| --- | --- | --- | --- | --- |") == 3)


    def row(i):
        r = [l for l in lines if l.startswith(f"| {i} |")]
        return r[0] if len(r) == 1 else None


    r1, r2, r3, r4 = row("z-one"), row("a-two"), row("z-three"), row("n-four")
    check("rows exist once each", all([r1, r2, r3, r4]), str([bool(x) for x in (r1, r2, r3, r4)]))
    if all([r1, r2, r3, r4]):
        c = cells(r1)
        check("row z-one has five cells", len(c) == 5, str(c))
        if len(c) == 5:
            check("z-one links", "[page](https://zeta.example/reports)" in c[2] and "[file](https://zeta.example/files/one.pdf)" in c[2], c[2])
            check("z-one save as", c[3] == "`ONE.pdf`", c[3])
            check("z-one notes first sentence only", c[4] == "First sentence here.", c[4])
            check("z-one title", c[1] == "Report One", c[1])
        c = cells(r2)
        check("row a-two has five cells (the pipe is escaped)", len(c) == 5, str(c))
        if len(c) == 5:
            check("a-two title escaped and unicode kept", c[1] == "Pipe \\| title — café", c[1])
            check("a-two only a page link", c[2] == "[page](https://alpha.example/)", c[2])
            check("a-two warning shown first", c[4].startswith("warning: filename_year_mismatch"), c[4])
            check("a-two note is one line", "\n" not in c[4] and "Line one line two." in c[4], c[4])
        c = cells(r3)
        if len(c) == 5:
            check("z-three link only and no notes", c[3] == "link only" and c[4] == "", str(c))
        c = cells(r4)
        if len(c) == 5:
            check("n-four long note truncated to 160 characters with an ellipsis", len(c[4]) <= 160 and c[4].endswith("…") and "Tail sentence" not in c[4], f"{len(c[4])} {c[4][-30:]!r}")
        check("catalog order kept inside a publisher", lines.index(r1) < lines.index(r3))
        check("publisher order in the file", lines.index("### alpha Foundation") < lines.index(r2) < lines.index("### Unknown publisher") < lines.index(r4) < lines.index("### Zeta Org") < lines.index(r1))

    check("request line", "- **Ask Registry** (Registry Hospital): [page](https://registry.example/) - needs permission or an application; do not automate." in lines)
    check("forbidden line", "- **Closed Database** (Closed DB): [page](https://closed.example/) - the publisher refuses automated access." in lines)
    for absent in ("fetch-x", "Fetchable report", "api-y", "An API", "X.pdf"):
        check(f"fetch and api entries are not listed ({absent})", absent not in md)
    check("prompt block exact", f"```text\n{PROMPT}\n```" in md)
    check("prompt is last", md.rstrip("\n").endswith("```"))

    empty = mod.render({"schema_version": 1, "sources": [{"id": "x", "access": "fetch", "title": "t"}]})
    check("empty catalog count line", "0 to download by hand, 0 to ask for, 0 not to download." in empty.split("\n"))
    check("empty catalog has no browser section", "## Download in a browser" not in empty and "## Ask first (access by request)" not in empty and "## Do not download" not in empty)
    check("empty catalog still has the prompt", PROMPT in empty)

    # command line
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        cat = tmp / "catalog.yaml"
        cat.write_text(yaml.safe_dump(CAT, allow_unicode=True, sort_keys=False), encoding="utf-8")
        out = tmp / "out" / "manual-downloads.md"
        out.parent.mkdir()

        def run(*args):
            return subprocess.run([sys.executable, str(cand), "--catalog", str(cat), "--out", str(out), *args], capture_output=True, text=True, encoding="utf-8", stdin=subprocess.DEVNULL, timeout=60)

        p = run("--check")
        check("--check on a missing file exits 1", p.returncode == 1 and "out of date" in p.stdout, f"{p.returncode} {p.stdout!r} {p.stderr[:100]!r}")
        p = run()
        check("write exits 0", p.returncode == 0, f"{p.returncode} {p.stderr[:160]!r}")
        if out.exists():
            data = out.read_bytes()
            check("file bytes equal render() as UTF-8 with LF", data == md.encode("utf-8"))
            check("no BOM", not data.startswith(b"\xef\xbb\xbf"))
        else:
            check("file written", False)
        p = run("--check")
        check("--check on a fresh file exits 0", p.returncode == 0 and "up to date" in p.stdout, f"{p.returncode} {p.stdout!r}")
        if out.exists():
            out.write_bytes(out.read_bytes() + b"x")
        p = run("--check")
        check("--check on a changed file exits 1 and writes nothing", p.returncode == 1 and out.exists() and out.read_bytes().endswith(b"x"), f"{p.returncode}")
        p = run()
        check("running again repairs the file", p.returncode == 0 and out.exists() and out.read_bytes() == md.encode("utf-8"))
        # Git on Windows (core.autocrlf) checks the page out with CRLF while the committed bytes are LF: --check must accept that, and still reject a stale page
        out.write_bytes(md.encode("utf-8").replace(b"\n", b"\r\n"))
        p = run("--check")
        check("--check accepts the same page checked out with CRLF line endings", p.returncode == 0 and "up to date" in p.stdout, f"{p.returncode} {p.stdout!r}")
        out.write_bytes(md.encode("utf-8").replace(b"\n", b"\r\n") + b"stale line\r\n")
        p = run("--check")
        check("--check still rejects a stale page that has CRLF line endings", p.returncode == 1 and "out of date" in p.stdout, f"{p.returncode} {p.stdout!r}")
        p = run()
        check("after the CRLF checks, running again rewrites the LF page", p.returncode == 0 and out.read_bytes() == md.encode("utf-8"))

    # the real catalog and the page committed with it
    REPO = Path(__file__).resolve().parents[2]
    real = REPO / "sources" / "catalog.yaml"
    page = REPO / "sources" / "manual-downloads.md"
    if real.exists() and page.exists():
        rc = yaml.safe_load(real.read_text(encoding="utf-8"))
        rendered = mod.render(rc)
        # Git on Windows (core.autocrlf) checks the page out with CRLF; the committed bytes are LF.
        check("the committed page is current (run make_manual_downloads.py)", page.read_bytes().replace(b"\r\n", b"\n") == rendered.encode("utf-8"))
        for s in rc["sources"]:
            listed = s.get("title", s["id"]) in rendered or s["id"] in rendered
            if s["access"] in ("manual", "request", "forbidden"):
                check(f"{s['id']} is on the page", listed)
            else:
                check(f"{s['id']} ({s['access']}) is not on the page", not listed)
    return fails


class Checks(unittest.TestCase):
    def test_every_check_passes(self):
        found = run_checks()
        self.assertEqual(found, [], "\n".join(found))


if __name__ == "__main__":
    # Run directly, an optional first argument names a candidate script to test instead of the real one.
    fails = run_checks(sys.argv[1] if len(sys.argv) > 1 else None)
    print("\n".join(fails[:25]))
    print("VERIFIED" if not fails else f"{len(fails)} failing check(s)")
    sys.exit(1 if fails else 0)
