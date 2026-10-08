"""Test for check_catalog_fields.py (board T-0026): runs it on small fixtures. Usage: python test_check_catalog_fields.py [SCRIPT]
Every failure line shows what was run, what was expected and what came back. Prints VERIFIED and exits 0 when all pass."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

DEFAULT_SCRIPT = Path(__file__).resolve().parent / "check_catalog_fields.py"

CATALOG_OK = """schema_version: 1
updated: 2026-10-04
defaults:
  redistribute: false
sources:
  - id: a
    title: A
    publisher: P
    kind: k
    access: fetch
    url: https://example.invalid/a
    notes: n
  - id: b
    title: B
    publisher: P
    kind: k
    access: api
    base_url: https://example.invalid/
    terms_url: https://example.invalid/t
"""
README_OK = """# Sources

Some text with a table that is NOT the fields table:

| `access` | Meaning |
| --- | --- |
| `fetch` | download |
| `api` | call it |

## Fields

| Field | Meaning |
| --- | --- |
| `id` | slug |
| `title`, `publisher`, `kind` | what it is |
| `access` | permission |
| `url` | direct url |
| `notes` | free text |
| `base_url` | api root |
| `terms_url` | usage policy |
| `redistribute` | set once under defaults |

## After the table
| `bogus` | not part of the fields table |
"""


def run_checks(script=None):
    """Run every check against SCRIPT (default: the real script next to this file). Returns the FAIL lines; empty when all pass."""
    cand = Path(script).resolve() if script else DEFAULT_SCRIPT
    fails = []



    def run(catalog, readme, args=()):
        d = Path(tempfile.mkdtemp(prefix="cf "))
        c, r = d / "catalog.yaml", d / "README.md"
        if catalog is not None:
            c.write_text(catalog, encoding="utf-8")
        if readme is not None:
            r.write_text(readme, encoding="utf-8")
        p = subprocess.run([sys.executable, str(cand), "--catalog", str(c), "--readme", str(r), *args], capture_output=True, text=True, encoding="utf-8", timeout=60, stdin=subprocess.DEVNULL)
        return p.returncode, (p.stdout + p.stderr)


    def expect(name, got, want_code, must_contain=(), must_not_contain=()):
        code, out = got
        if code != want_code:
            fails.append(f"FAIL {name}: expected exit {want_code}, got {code}. Output was: {out.strip()[-300:]!r}")
            return
        for s in must_contain:
            if s not in out:
                fails.append(f"FAIL {name}: expected the output to contain {s!r}. Output was: {out.strip()[-300:]!r}")
        for s in must_not_contain:
            if s in out:
                fails.append(f"FAIL {name}: the output must not contain {s!r}. Output was: {out.strip()[-300:]!r}")


    try:
        expect("all documented", run(CATALOG_OK, README_OK), 0, ["fields ok"])
        expect("a key used in the catalog but missing from the table", run(CATALOG_OK.replace("notes: n", "notes: n\n    brand_new: 1"), README_OK), 1, ["undocumented: brand_new"])
        expect("two missing keys are both named", run(CATALOG_OK.replace("notes: n", "notes: n\n    zeta: 1\n    alpha: 2"), README_OK), 1, ["undocumented: alpha, zeta"])
        expect("a key listed in the table but never used is a warning only", run(CATALOG_OK, README_OK.replace("| `notes` | free text |", "| `notes` | free text |\n| `legacy_field` | old |")), 0, ["unused (warning): legacy_field", "fields ok"])
        expect("a defaults key counts as documented when listed", run(CATALOG_OK, README_OK), 0, [], ["undocumented"])
        expect("a defaults key missing from the table is reported", run(CATALOG_OK, README_OK.replace("| `redistribute` | set once under defaults |\n", "")), 1, ["undocumented: redistribute"])
        expect("the access-values table is not the fields table", run(CATALOG_OK, README_OK), 0, [], ["unused (warning): fetch", "unused (warning): api"])
        expect("rows after the fields table are ignored", run(CATALOG_OK, README_OK), 0, [], ["bogus"])
        expect("a missing README exits 2", run(CATALOG_OK, None), 2, ["error:"])
        expect("a README without a fields table exits 2", run(CATALOG_OK, "# nothing here\n"), 2, ["error:"])
        expect("an invalid catalog exits 2", run("sources: [unclosed", README_OK), 2, ["error:"])
        expect("a catalog without a sources list exits 2", run("schema_version: 1\n", README_OK), 2, ["error:"])
        code, out = run(CATALOG_OK, README_OK, ["--bogus-flag"])
        if code == 0:
            fails.append("FAIL an unknown flag must not be accepted silently (argparse exits 2)")
    except Exception as e:  # noqa: BLE001
        fails.append(f"FAIL the verifier could not finish: {type(e).__name__}: {e}")
    return fails


class Checks(unittest.TestCase):
    def test_every_check_passes(self):
        found = run_checks()
        self.assertEqual(found, [], "\n".join(found))


if __name__ == "__main__":
    # Run directly, an optional first argument names a candidate script to test instead of the real one.
    fails = run_checks(sys.argv[1] if len(sys.argv) > 1 else None)
    for f in fails:
        print(f)
    print("VERIFIED" if not fails else f"{len(fails)} failure(s)")
    sys.exit(1 if fails else 0)
