"""Re-stamp content_hash and confirmed_on (and, with --confirmed-by, confirmed_by) on rows whose hash no longer matches
only because summary_es was added or changed. Line-preserving; touches only those lines per row.

Usage: python reconfirm_summary_es.py --yaml <cf-research>/resources/resources.yaml --date YYYY-MM-DD
           [--ids id1,id2] [--confirmed-by council-on-behalf-of-maintainer] [--base-yaml OLD.yaml]

Without --confirmed-by, confirmed_by is left as it is: that is the maintainer's own path, run AFTER reading every
summary_es. --confirmed-by accepts only council-on-behalf-of-maintainer (the Spanish was approved by a model council
under the maintainer's written authorisation); this script can never write confirmed_by: maintainer.

Refuses a row whose hash does not match its OLD content (name, url, summary, who_for without summary_es), so it can
never bless a change to any field other than summary_es. A row whose already-confirmed summary_es was changed is
accepted only when it is named in --ids and --base-yaml holds the previous version of the file: the previous row must
match its own content_hash, and name, url, summary and who_for must be unchanged since then.
"""
import argparse, hashlib, re, sys
import yaml

# Values this script may write to confirmed_by. "maintainer" is deliberately absent: only a person writes it.
SCRIPT_CONFIRMED_BY = ("council-on-behalf-of-maintainer",)
BASE_FIELDS = ("name", "url", "summary", "who_for")


def h(row, with_es):
    parts = [str(row.get(k, "")) for k in BASE_FIELDS]
    if with_es and "summary_es" in row:
        parts.append(str(row["summary_es"]))
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--date", required=True)
    ap.add_argument("--ids", default="")
    ap.add_argument("--confirmed-by", default=None, choices=SCRIPT_CONFIRMED_BY,
                    help="also set confirmed_by (default: leave it unchanged)")
    ap.add_argument("--base-yaml", default=None,
                    help="previous version of the file, to accept explicitly named rows whose summary_es changed")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    raw = open(a.yaml, "rb").read().decode("utf-8")
    rows = {r["id"]: r for r in yaml.safe_load(raw)["resources"]}
    named = set(filter(None, a.ids.split(",")))
    base = {}
    if a.base_yaml:
        base = {r["id"]: r for r in yaml.safe_load(open(a.base_yaml, encoding="utf-8"))["resources"]}
    want = named or {
        i for i, r in rows.items() if "summary_es" in r and "content_hash" in r and r["content_hash"] != h(r, True)}
    for i in sorted(want):
        if i not in rows:
            sys.exit(f"refused: {i} is not in the file")
        r = rows[i]
        if "summary_es" not in r or "content_hash" not in r:
            sys.exit(f"refused: {i} has no summary_es or no content_hash")
        if r["content_hash"] == h(r, False):
            continue
        old = base.get(i)
        if (i in named and old is not None and old.get("content_hash") == h(old, True)
                and r["content_hash"] == old["content_hash"]
                and all(str(old.get(k, "")) == str(r.get(k, "")) for k in BASE_FIELDS)):
            continue
        sys.exit(f"refused: {i} changed in more than summary_es; confirm it by hand")
    out, cur = [], None
    for ln in raw.splitlines(keepends=True):
        m = re.match(r"^  - id: (\S+)\s*$", ln)
        if m:
            cur = m.group(1)
        if cur in want:
            eol = "\r\n" if ln.endswith("\r\n") else "\n"
            if ln.startswith("    content_hash: "):
                ln = f"    content_hash: {h(rows[cur], True)}{eol}"
            elif ln.startswith("    confirmed_on: "):
                ln = f"    confirmed_on: {a.date}{eol}"
            elif a.confirmed_by and ln.startswith("    confirmed_by: "):
                ln = f"    confirmed_by: {a.confirmed_by}{eol}"
        out.append(ln)
    open(a.yaml, "wb").write("".join(out).encode("utf-8"))
    who = f", confirmed_by={a.confirmed_by}" if a.confirmed_by else ""
    print(f"re-stamped {len(want)} rows (content_hash, confirmed_on={a.date}{who})")


if __name__ == "__main__":
    main()
