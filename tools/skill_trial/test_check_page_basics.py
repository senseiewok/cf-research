"""Test for check_page_basics.py: a checker that has never failed proves nothing.

A good page must pass every check; a good page with a properly guarded animation must too; and each bad stub is the good page with ONE change
(made here by text replacement) and must fail the check that change breaks. Needs a Chromium or Edge; no network.
Usage: python test_check_page_basics.py     One FAIL line per failing check; VERIFIED and exit 0 when all pass."""
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKER = HERE / "check_page_basics.py"
fails = []

COPY = [
    "Sensei Ewok's CF Lab builds reusable AI and MCP tools in support of cystic fibrosis research, alongside creative computing.",
    "The work is open source and done in donated time.",
    "Nothing here is medical advice or a claim about any treatment.",
]

GOOD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>What this lab is</title>
<style>
  body { margin: 0; background: #fafaf7; color: #1d2430; font: 1rem/1.6 Georgia, serif; }
  main { max-width: 38rem; margin: 0 auto; padding: 1.5rem 1rem; }
  h1 { font-size: 1.6rem; line-height: 1.2; }
  a { color: #0b4f8a; }
</style>
</head>
<body>
<main>
  <h1>What this lab is</h1>
  <p>Sensei Ewok's CF Lab builds reusable AI and MCP tools in support of cystic fibrosis research, alongside creative computing.</p>
  <p>The work is open source and done in donated time.</p>
  <p>Nothing here is medical advice or a claim about any treatment.</p>
</main>
</body>
</html>
"""

GOOD_DATA_URI = GOOD.replace("</style>", "  body::before { content: ''; display: block; height: 4px; background: url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence/%3E%3C/filter%3E%3Crect filter='url(%23n)' width='4' height='4'/%3E%3C/svg%3E\"); }\n</style>")

GOOD_ANIMATED = GOOD.replace("</style>", "  .fade { animation: fade 0.8s ease forwards; }\n  @keyframes fade { from { opacity: 0.2; } to { opacity: 1; } }\n  @media (prefers-reduced-motion: reduce) { .fade { animation: none; } }\n</style>").replace("<h1>", "<h1 class=\"fade\">")

# name: (the change as (old, new) pairs, the check it must fail)
STUBS = {
    "no lang attribute": ([('<html lang="en">', "<html>")], "lang"),
    "two h1 headings": ([("</main>", "<h1>Another</h1></main>")], "lang"),
    "no main landmark": ([("<main>", "<div>"), ("</main>", "</div>")], "lang"),
    "empty title": ([("<title>What this lab is</title>", "<title></title>")], "lang"),
    "low contrast text": ([("color: #1d2430;", "color: #8a8a82;")], "contrast"),
    "text over a gradient": ([("background: #fafaf7;", "background: linear-gradient(#ffffff, #eeeeee);")], "contrast"),
    "scrolls sideways at 360 px": ([("</main>", '<div style="width:600px;height:8px;background:#000"></div></main>')], "mobile"),
    "an external stylesheet": ([("<style>", '<link rel="stylesheet" href="https://example.com/x.css"><style>')], "deps"),
    "an external script": ([("</body>", '<script src="https://example.com/x.js"></script></body>')], "deps"),
    "an @import of a font": ([("<style>", '<style>@import url("https://example.com/f.css");')], "deps"),
    "a css background image from elsewhere": ([("</style>", "  main { background: url(https://example.com/bg.png); }\n</style>")], "deps"),
    "an image from elsewhere": ([("</main>", '<img src="https://example.com/a.png" alt="a"></main>')], "deps"),
    "an animation that runs forever": ([("</style>", "  .p { animation: pulse 1s infinite; }\n  @keyframes pulse { from { opacity: 0.5; } to { opacity: 1; } }\n  @media (prefers-reduced-motion: reduce) { .p { animation: none; } }\n</style>"), ("<h1>", '<h1 class="p">')], "motion"),
    "an animation with no reduced-motion rule": ([("</style>", "  .fade { animation: fade 0.8s ease forwards; }\n  @keyframes fade { from { opacity: 0.2; } to { opacity: 1; } }\n</style>"), ("<h1>", '<h1 class="fade">')], "motion"),
    "a verified sentence changed": ([("done in donated time", "done in spare time")], "copy"),
    "a verified sentence missing": ([("  <p>Nothing here is medical advice or a claim about any treatment.</p>\n", "")], "copy"),
    "a console error": ([("</body>", "<script>undefinedFunction();</script></body>")], "console"),
}


def check(name, ok, detail=""):
    if not ok:
        fails.append(name)
        print(f"FAIL {name}" + (f"  ({detail})" if detail else ""))


def run(page, copy):
    p = subprocess.run([sys.executable, str(CHECKER), str(page), "--copy", str(copy)], capture_output=True, text=True, timeout=300)
    return p.returncode, p.stdout, p.stderr


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    copy = tmp / "copy.txt"
    copy.write_text("\n".join(COPY) + "\n", encoding="utf-8")
    for label, src in (("the good page", GOOD), ("the good page with a guarded animation", GOOD_ANIMATED), ("the good page with an inline SVG image that holds a url(#...) reference", GOOD_DATA_URI)):
        page = tmp / (label.replace(" ", "-") + ".html")
        page.write_text(src, encoding="utf-8")
        code, out, err = run(page, copy)
        last = out.strip().splitlines()[-1] if out.strip() else ""
        check(f"{label} passes every check", code == 0 and last == "7/7 checks passed", f"exit={code} out={out[:400]!r} err={err[:120]!r}")
    for name, (edits, want) in STUBS.items():
        src = GOOD
        applied = True
        for old, new in edits:
            if old not in src:
                applied = False
                break
            src = src.replace(old, new, 1)
        if not applied:
            check(f"stub '{name}': the edit applies to the good page", False, "text to replace not found")
            continue
        page = tmp / f"stub-{abs(hash(name))}.html"
        page.write_text(src, encoding="utf-8")
        code, out, err = run(page, copy)
        failed = [l.split(":")[0].replace("FAIL ", "") for l in out.splitlines() if l.startswith("FAIL ")]
        check(f"stub '{name}': exits 1", code == 1, f"exit={code}")
        check(f"stub '{name}': fails the '{want}' check", want in failed, f"failed: {failed}; out={out[:300]!r}")
    code, out, err = run(tmp / "nope.html", copy)
    check("a missing page exits 2", code == 2 and out.startswith("ERROR"), f"exit={code} {out[:80]!r}")

print("VERIFIED" if not fails else f"{len(fails)} failing check(s)")
sys.exit(1 if fails else 0)
