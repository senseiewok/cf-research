"""For the maintainer, AFTER reading every summary_es: re-stamp content_hash (and confirmed_on) on rows whose hash no
longer matches only because summary_es was added. Line-preserving; touches only those two lines per row.

Usage: python reconfirm_summary_es.py --yaml <cf-research>/resources/resources.yaml --date YYYY-MM-DD [--ids id1,id2]
Refuses a row whose hash does not match its OLD content (name, url, summary, who_for without summary_es), so it can
never bless a change to any field other than summary_es.
"""
import argparse, hashlib, re, sys
import yaml


def h(row, with_es):
    parts = [str(row.get(k, "")) for k in ("name", "url", "summary", "who_for")]
    if with_es and "summary_es" in row:
        parts.append(str(row["summary_es"]))
    return hashlib.sha256("\n".join(parts).encode("utf-8")).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--date", required=True)
    ap.add_argument("--ids", default="")
    a = ap.parse_args()
    raw = open(a.yaml, "rb").read().decode("utf-8")
    rows = {r["id"]: r for r in yaml.safe_load(raw)["resources"]}
    want = set(filter(None, a.ids.split(","))) or {
        i for i, r in rows.items() if "summary_es" in r and "content_hash" in r and r["content_hash"] != h(r, True)}
    for i in sorted(want):
        r = rows[i]
        if r["content_hash"] != h(r, False):
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
        out.append(ln)
    open(a.yaml, "wb").write("".join(out).encode("utf-8"))
    print(f"re-stamped {len(want)} rows (content_hash, confirmed_on={a.date})")


if __name__ == "__main__":
    main()
