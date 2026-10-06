"""Test for check_numbers.py: every number in a draft must appear in the evidence, or be allowed with a written reason.
Usage: python test_check_numbers.py [SCRIPT]   One FAIL line per failing check; prints VERIFIED and exits 0 when all pass."""
import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path

cand = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent / "check_numbers.py"
try:
    spec = importlib.util.spec_from_file_location("cand_cn", cand)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
except Exception as e:  # noqa: BLE001
    print(f"FAIL import: {type(e).__name__}: {str(e)[:160]}")
    sys.exit(1)

fails = []


def check(name, cond, detail=""):
    if not cond:
        fails.append(f"FAIL {name}: {detail}")


EVIDENCE = "In 2024, 99.4 percent (n=33,782) were genotyped; 1,214 variants; F508del 85.0% and R117H 3.5 percent; 2,000 variants."


def bad(draft, evidence=EVIDENCE, allow=None, ignore_below=10):
    return [h["number"] for h in mod.unmatched_numbers(draft, evidence, allow or {}, ignore_below)]


check("numbers in the evidence pass", bad("99.4 percent of 33,782 people, 1,214 variants, 85.0% and 3.5 percent in 2024.") == [])
check("an invented number is flagged", bad("About 99.7 percent were genotyped.") == ["99.7"], str(bad("About 99.7 percent were genotyped.")))
check("a changed digit is flagged", bad("There were 33,785 people.") == ["33785"])
check("commas and trailing zeros do not matter", bad("33782 people and 99.40 percent and 85 percent and 86 percent") == ["86"], str(bad("33782 people and 99.40 percent and 85 percent and 86 percent")))
check("a percent sign is not part of the number", bad("85.0% of them") == [])
check("small integers are ignored by default", bad("There were 3 groups and 7 reports and 9 checks.") == [])
check("ten and above are checked", bad("There were 10 groups.") == ["10"])
check("decimals are always checked, even small", bad("A ratio of 4.4 was seen.") == ["4.4"])
check("ignore_below=0 checks small integers", bad("There were 3 groups.", ignore_below=0) == ["3"])
check("dates and ids are not numbers", bad("Read on 2026-10-05, board T-0067, variant R117H, build hg19, NC_000007.13, version 1.0.0, doi 10.1002/humu.23276.") == [])
check("units attached to a word are not numbers", bad("A 40s age band and 5T alleles and 2x faster.") == [])
check("code is ignored", bad("Run `check --max-run 25` now.\n```\nx = 98765\n```\n") == [])
check("each unmatched number is reported once per line with its line", [h["line"] for h in mod.unmatched_numbers("ok\n\nthe 99.7 case\n", EVIDENCE, {}, 10)] == [3])
check("the same number twice on a line is reported once", len(mod.unmatched_numbers("99.7 and again 99.7", EVIDENCE, {}, 10)) == 1)
check("the same number on two lines is reported on each", [h["line"] for h in mod.unmatched_numbers("99.7\n99.7", EVIDENCE, {}, 10)] == [1, 2])
check("allowed numbers pass", bad("A ratio of 4.4 and 564 in total.", allow={"4.4": "564 / 127", "564": "127 + 437"}) == [])
check("an allowed number must match exactly", bad("A ratio of 4.5.", allow={"4.4": "564 / 127"}) == ["4.5"])
check("empty draft and empty evidence", bad("") == [] and bad("There were 1,234 cases.", evidence="") == ["1234"])

with tempfile.TemporaryDirectory() as td:
    d = Path(td)
    (d / "ev.txt").write_text(EVIDENCE, encoding="utf-8")
    (d / "ok.md").write_text("99.4 percent were genotyped.", encoding="utf-8")
    (d / "bad.md").write_text("99.7 percent were genotyped.\nAnd 1,214 variants.", encoding="utf-8")
    (d / "allow.txt").write_text("99.7  # a deliberately allowed test number\n# a comment line\n\n", encoding="utf-8")
    (d / "noreason.txt").write_text("99.7\n", encoding="utf-8")
    (d / "emptyreason.txt").write_text("99.7  #   \n", encoding="utf-8")

    def run(*args):
        return subprocess.run([sys.executable, str(cand), *args], capture_output=True, text=True)

    r = run("--evidence", str(d / "ev.txt"), str(d / "ok.md"))
    check("cli exit 0", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    r = run("--evidence", str(d / "ev.txt"), str(d / "bad.md"))
    check("cli exit 1 names file line and number", r.returncode == 1 and "bad.md" in r.stdout and "99.7" in r.stdout and ":1:" in r.stdout, f"{r.returncode} {r.stdout[-300:]}")
    check("cli does not flag the number that is in the evidence", "1,214" not in r.stdout and "1214" not in r.stdout, r.stdout[-300:])
    r = run("--evidence", str(d / "ev.txt"), "--allow-file", str(d / "allow.txt"), str(d / "bad.md"))
    check("cli allow-file with a reason passes", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    for f in ("noreason.txt", "emptyreason.txt"):
        r = run("--evidence", str(d / "ev.txt"), "--allow-file", str(d / f), str(d / "bad.md"))
        check(f"cli allow-file without a reason is an error ({f})", r.returncode == 2 and "reason" in r.stdout, f"{r.returncode} {r.stdout[-200:]}")
    r = run("--evidence", str(d / "nope.txt"), str(d / "ok.md"))
    check("cli missing evidence exits 2", r.returncode == 2, f"{r.returncode}")
    r = run("--evidence", str(d / "ev.txt"), str(d / "nope.md"))
    check("cli missing draft exits 2", r.returncode == 2, f"{r.returncode}")
    r = run(str(d / "ok.md"))
    check("cli with no evidence at all exits 2", r.returncode == 2, f"{r.returncode} {r.stdout[-200:]}")
    r = run("--evidence", str(d), str(d / "ok.md"))
    check("cli reads a folder of text files as evidence", r.returncode == 0, f"{r.returncode} {r.stdout[-200:]}")
    (d / "u.md").write_text("The count ✓ was 99.7 ☃.", encoding="utf-8")
    r = subprocess.run([sys.executable, str(cand), "--evidence", str(d / "ev.txt"), str(d / "u.md")], capture_output=True, env=dict(os.environ, PYTHONIOENCODING="cp1252"))
    check("cli survives characters the console cannot print", r.returncode == 1 and b"Traceback" not in r.stderr, f"{r.returncode} {r.stderr[-200:]}")

if fails:
    print("\n".join(fails))
    sys.exit(1)
print("VERIFIED")
