"""Wrap the lab's check-ascii.py so a trial can judge ASCII-art output without being fooled by a vacuous pass.

check-ascii.py exits 0 for a file with NO art block at all ("0 art block(s) checked: 0 error(s)"). A trial that counted that as a pass would reward an
answer that draws nothing, so this wrapper requires at least one block, and fails on warnings too (--strict). It prints one PASS or FAIL line per check:
  block   at least one art block was found and checked
  rules   check-ascii.py reports no error and no warning at the given width

Usage: python check_ascii_art.py FILE.md --ascii-checker PATH/check-ascii.py [--max-width 36]
Exit 0 when both pass, 1 otherwise, 2 for a missing file or checker."""
import argparse
import re
import subprocess
import sys
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("file")
    ap.add_argument("--ascii-checker", required=True)
    ap.add_argument("--max-width", type=int, default=36)
    a = ap.parse_args(argv)
    if not Path(a.file).is_file() or not Path(a.ascii_checker).is_file():
        print(f"ERROR: missing file or checker: {a.file} / {a.ascii_checker}")
        return 2
    p = subprocess.run([sys.executable, a.ascii_checker, "--strict", "--max-width", str(a.max_width), a.file], capture_output=True, text=True, timeout=120)
    out = p.stdout + p.stderr
    m = re.search(r"(\d+) art block\(s\) checked", out)
    blocks = int(m.group(1)) if m else 0
    findings = [ln.strip() for ln in out.splitlines() if re.search(r"\b(ERROR|WARN)\b", ln)]
    results = [("block", blocks >= 1, f"{blocks} art block(s) found; a file with none passes check-ascii.py vacuously"),
               ("rules", p.returncode == 0 and blocks >= 1, "; ".join(findings[:3]) if findings else f"exit {p.returncode}")]
    for name, ok, detail in results:
        print(("PASS " if ok else "FAIL ") + name + ("" if ok else f": {detail}"))
    passed = sum(1 for _n, ok, _d in results if ok)
    print(f"{passed}/{len(results)} checks passed")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
