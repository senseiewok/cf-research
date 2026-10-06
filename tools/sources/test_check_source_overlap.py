"""Test for check_source_overlap.py: a draft may share only short runs of words with the source documents it was written from.
Usage: python test_check_source_overlap.py [SCRIPT]
One FAIL line per failing check; prints VERIFIED and exits 0 when all pass."""
import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

cand = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent / "check_source_overlap.py"
try:
    spec = importlib.util.spec_from_file_location("cand_so", cand)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
except Exception as e:  # noqa: BLE001
    print(f"FAIL import: {type(e).__name__}: {str(e)[:160]}")
    sys.exit(1)

fails = []


def check(name, cond, detail=""):
    if not cond:
        fails.append(f"FAIL {name}: {detail}")


def words(n, tag="w"):
    return " ".join(f"{tag}{i}" for i in range(n))


SRC_A = "Intro text here. " + words(40, "alpha") + " closing words of the source."
SRC_B = "Another source entirely. " + words(40, "beta") + " end."

with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    src = d / "src"
    src.mkdir()
    (src / "AAA_report.txt").write_text(SRC_A, encoding="utf-8")
    (src / "UKCFR_report.txt").write_text(SRC_B, encoding="utf-8")
    (src / "notes.md").write_text(" ".join(f"gamma{i}" for i in range(14)), encoding="utf-8")

    def run(text, **kw):
        p = d / "draft.md"
        p.write_text(text, encoding="utf-8")
        return mod.check_document(p, src, **kw)

    # a run just under the limit passes; at the limit plus one it fails
    ok = run("Our own words. " + " ".join(f"alpha{i}" for i in range(25)) + " then more of our words.", max_run=25)
    check("run of 25 passes", not ok, str(ok))
    bad = run("Lead in. " + " ".join(f"alpha{i}" for i in range(26)) + " tail.", max_run=25)
    check("run of 26 fails", len(bad) == 1 and bad[0]["source"] == "AAA_report.txt" and bad[0]["words"] == 26, str(bad))

    # case, punctuation and line breaks do not hide a copy
    disguised = "\n".join(f"Alpha{i}," for i in range(30))
    check("punctuation and case ignored", len(run(disguised, max_run=25)) == 1)

    # a paraphrase passes; a copy broken every ten words by our own words passes at limit 25
    broken = " ".join(f"alpha{i}" if i % 10 else "ours" for i in range(40))
    check("copy interrupted by own words passes", not run(broken, max_run=25))

    # the stricter limit for a named source applies only to that source
    uk = " ".join(f"beta{i}" for i in range(15))
    check("default limit lets 15 words of UK text through", not run(uk, max_run=25))
    check("per-source limit catches 15 words", len(run(uk, max_run=25, limits={"UKCFR": 10})) == 1)
    check("per-source limit does not touch other sources", not run(" ".join(f"alpha{i}" for i in range(15)), max_run=25, limits={"UKCFR": 10}))

    # only text and PDF sources are read; the .md file in the source folder is ignored
    check("md source ignored", run(" ".join(f"gamma{i}" for i in range(14)), max_run=12) == [])

    # a run split across two sources is not one run
    check("two sources not joined", not run(" ".join(f"alpha{i}" for i in range(20)) + " " + " ".join(f"beta{i}" for i in range(20)), max_run=25))

    # an empty draft and a draft shorter than the limit
    check("empty draft", run("", max_run=25) == [])
    check("short draft", run("alpha0 alpha1", max_run=25) == [])

    # a missing source folder or an empty one is an error, never a silent pass
    try:
        mod.check_document(d / "draft.md", d / "nope", max_run=25)
        check("missing source dir raises", False, "no exception")
    except Exception as e:  # noqa: BLE001
        check("missing source dir raises", isinstance(e, (FileNotFoundError, ValueError)), type(e).__name__)
    empty = d / "empty"
    empty.mkdir()
    try:
        mod.check_document(d / "draft.md", empty, max_run=25)
        check("empty source dir raises", False, "no exception")
    except Exception as e:  # noqa: BLE001
        check("empty source dir raises", isinstance(e, ValueError), type(e).__name__)

    # command line: exit 1 with the file, source and snippet on a copy; exit 0 on clean text; exit 2 on a bad folder
    (d / "copy.md").write_text("x " + " ".join(f"alpha{i}" for i in range(30)), encoding="utf-8")
    (d / "clean.md").write_text("Our own sentence about registries.", encoding="utf-8")
    r = subprocess.run([sys.executable, str(cand), "--sources", str(src), str(d / "copy.md")], capture_output=True, text=True)
    check("cli exit 1 on copy", r.returncode == 1, f"{r.returncode} {r.stdout[-200:]} {r.stderr[-200:]}")
    check("cli names file and source", "copy.md" in r.stdout and "AAA_report.txt" in r.stdout, r.stdout[-300:])
    r = subprocess.run([sys.executable, str(cand), "--sources", str(src), str(d / "clean.md")], capture_output=True, text=True)
    check("cli exit 0 on clean", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    r = subprocess.run([sys.executable, str(cand), "--sources", str(d / "nope"), str(d / "clean.md")], capture_output=True, text=True)
    check("cli exit 2 on bad folder", r.returncode == 2, f"{r.returncode}")

    # the snippet does not reprint the copied run in full (a report must not itself republish the text)
    r = subprocess.run([sys.executable, str(cand), "--sources", str(src), str(d / "copy.md")], capture_output=True, text=True)
    check("report does not reprint the run", "alpha20" not in r.stdout, r.stdout[-300:])

# --- quotation length: no sentence-length quotations in a public note, whatever the source -------------------------------------
def q(text, n=8):
    return mod.find_quotations(text, n)


check("quote of 8 words passes", q('He said "one two three four five six seven eight" here.') == [])
long9 = q('He said "one two three four five six seven eight nine" here.')
check("quote of 9 words flagged", len(long9) == 1 and long9[0]["words"] == 9, str(long9))
check("curly quotes flagged", len(q("He said “one two three four five six seven eight nine ten” here.")) == 1)
check("quote wrapped over lines flagged", len(q('He said "one two three four\nfive six seven eight\nnine ten" here.')) == 1)
check("short labels pass", q('The class "unclassified" and the label "likely CFTR-RD" and "Not done".') == [])
check("two quotes pair correctly", len(q('A "x y" then "one two three four five six seven eight nine" end.')) == 1)
check("fenced code skipped", q('```\nprint("one two three four five six seven eight nine ten")\n```\n') == [])
check("inline code skipped", q('Run `say "one two three four five six seven eight nine ten"` now.') == [])
check("empty and no quotes", q("") == [] and q("plain words only") == [])
check("quote does not pair across paragraphs", q('An odd " one two three four.\n\nFive six seven eight nine ten " end.') == [])
check("snippet is only the start", "nine" not in long9[0]["snippet"], long9[0]["snippet"])
check("max 0 disables", q('"one two three four five six seven eight nine ten"', 0) == [])

with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    (d / "n.md").write_text('He said "one two three four five six seven eight nine ten" today.', encoding="utf-8")
    (d / "ok.md").write_text("No quotations at all.", encoding="utf-8")
    r = subprocess.run([sys.executable, str(cand), "--quotes-only", str(d / "n.md")], capture_output=True, text=True)
    check("cli quotes-only exit 1", r.returncode == 1 and "quotation" in r.stdout and "n.md" in r.stdout, f"{r.returncode} {r.stdout[-200:]}")
    check("cli quotes-only does not print the quotation", "nine" not in r.stdout, r.stdout[-200:])
    r = subprocess.run([sys.executable, str(cand), "--quotes-only", str(d / "ok.md")], capture_output=True, text=True)
    check("cli quotes-only exit 0", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    r = subprocess.run([sys.executable, str(cand), "--quotes-only", "--sources", str(d / "nope"), str(d / "ok.md")], capture_output=True, text=True)
    check("cli quotes-only needs no sources folder", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    r = subprocess.run([sys.executable, str(cand), "--quotes-only", "--max-quote", "20", str(d / "n.md")], capture_output=True, text=True)
    check("cli max-quote raises the limit", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    r = subprocess.run([sys.executable, str(cand), "--sources", str(d / "nope"), str(d / "n.md")], capture_output=True, text=True)
    check("cli with sources still exits 2 on a bad folder", r.returncode == 2, f"{r.returncode}")
    # a snippet with characters the console cannot show must not crash the report (Windows consoles use cp1252)
    (d / "u.md").write_text('He said "café ✓ ☃ two three four five six seven eight nine ten" today.', encoding="utf-8")
    env = dict(os.environ, PYTHONIOENCODING="cp1252")
    r = subprocess.run([sys.executable, str(cand), "--quotes-only", str(d / "u.md")], capture_output=True, env=env)
    check("cli survives characters the console cannot print", r.returncode == 1 and b"Traceback" not in r.stderr, f"{r.returncode} {r.stderr[-200:]}")

if fails:
    print("\n".join(fails))
    sys.exit(1)
print("VERIFIED")
