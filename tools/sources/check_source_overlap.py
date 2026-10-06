"""Check that a draft shares only short runs of words with the source documents it was written from.

The catalog records what each publisher permits (sources/catalog.yaml); most reports may be read, cited and linked but not reproduced, and
the UK report says its content may not be reproduced in publications without permission. This tool measures one thing a script can: the
longest run of consecutive words a draft has in common with any one source file. It does not decide what is fair quotation; a person does.

A second check needs no sources: a quoted span of more than --max-quote words (default 8) in a draft is a quotation of someone's sentence,
and the lab's rule is to paraphrase; names, titles, labels and field names are short enough to pass. `--quotes-only` runs just that check.

Usage: python check_source_overlap.py [--sources DIR] [--max-run N] [--limit PREFIX=N ...] [--max-quote N] [--quotes-only] DRAFT [DRAFT ...]
Sources are the .txt and .pdf files in DIR (default: sources/downloads, which is git-ignored); PDFs are read with pdftotext. Case, punctuation
and line breaks are ignored. A run is counted only where the words line up in the same order in one source. The report names the draft, the
source, the length and the first words of the run, never the whole run. Exit 0: no run over the limit. 1: at least one. 2: bad input."""
import argparse
import re
import subprocess
import sys
from pathlib import Path

DEFAULT_MAX_RUN = 25
# Stricter limits for sources whose own text asks for permission before reproduction in publications (catalog: ukcfr-adr-2024).
DEFAULT_LIMITS = {"UKCFR": 10}
SNIPPET_WORDS = 6
DEFAULT_MAX_QUOTE = 8
_FENCE = re.compile(r"^```.*?^```", re.S | re.M)
_INLINE_CODE = re.compile(r"`[^`\n]*`")
_QUOTED = re.compile(r'"([^"]*)"|“([^”]*)”')


def tokens(text):
    return re.findall(r"\w+", text.lower())


def source_text(path):
    if path.suffix.lower() == ".txt":
        return path.read_text(encoding="utf-8", errors="replace")
    try:
        r = subprocess.run(["pdftotext", str(path), "-"], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180)
    except FileNotFoundError as e:
        raise ValueError("pdftotext is not installed; give a folder of .txt files or install poppler") from e
    if r.returncode != 0:
        raise ValueError(f"pdftotext failed on {path.name}")
    return r.stdout


def limit_for(name, max_run, limits):
    for prefix, n in (limits or {}).items():
        if name.startswith(prefix):
            return n
    return max_run


def runs_over(doc, src, limit):
    """Runs of more than `limit` words in `doc` that appear in the same order in `src`. Returns [(start, length)]."""
    n = limit + 1
    if len(doc) < n or len(src) < n:
        return []
    index = {}
    for p in range(len(src) - n + 1):
        index.setdefault(tuple(src[p:p + n]), []).append(p)
    found, active, start = [], set(), 0
    for i in range(len(doc) - n + 1):
        here = {p - i for p in index.get(tuple(doc[i:i + n]), ())}
        keep = active & here
        if keep:
            active = keep
            continue
        if active:
            found.append((start, i - start + n - 1))
        active, start = here, i
    if active:
        found.append((start, len(doc) - n + 1 - start + n - 1))
    return found


def check_document(doc_path, sources_dir, max_run=DEFAULT_MAX_RUN, limits=None):
    sources_dir = Path(sources_dir)
    if not sources_dir.is_dir():
        raise FileNotFoundError(f"no such folder: {sources_dir}")
    files = sorted(p for p in sources_dir.iterdir() if p.is_file() and p.suffix.lower() in (".txt", ".pdf"))
    if not files:
        raise ValueError(f"no .txt or .pdf files in {sources_dir}")
    doc = tokens(Path(doc_path).read_text(encoding="utf-8", errors="replace"))
    out = []
    for f in files:
        limit = limit_for(f.name, max_run, limits)
        if len(doc) <= limit:
            continue
        for start, length in runs_over(doc, tokens(source_text(f)), limit):
            out.append({"source": f.name, "words": length, "limit": limit, "doc_word_index": start,
                        "snippet": " ".join(doc[start:start + SNIPPET_WORDS])})
    return out


def find_quotations(text, max_words=DEFAULT_MAX_QUOTE):
    """Quoted spans of more than `max_words` words, as [{"words", "snippet"}]. Code fences and inline code are ignored, and a quotation mark
    pairs only within its own paragraph, so one stray mark cannot swallow the text after it. 0 or less disables the check."""
    if max_words <= 0:
        return []
    text = _INLINE_CODE.sub("", _FENCE.sub("", text))
    found = []
    for para in re.split(r"\n\s*\n", text):
        for m in _QUOTED.finditer(para):
            words = (m.group(1) if m.group(1) is not None else m.group(2)).split()
            if len(words) > max_words:
                found.append({"words": len(words), "snippet": " ".join(words[:SNIPPET_WORDS])})
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description="Longest run of words a draft shares with a source document.")
    ap.add_argument("drafts", nargs="+")
    ap.add_argument("--max-quote", type=int, default=DEFAULT_MAX_QUOTE, help="flag a quoted span of more than N words (0 turns the check off)")
    ap.add_argument("--quotes-only", action="store_true", help="check quotation length only; needs no sources folder")
    ap.add_argument("--sources", default=str(Path(__file__).resolve().parents[2] / "sources" / "downloads"))
    ap.add_argument("--max-run", type=int, default=DEFAULT_MAX_RUN)
    ap.add_argument("--limit", action="append", default=[], metavar="PREFIX=N", help="stricter limit for source files whose name starts with PREFIX")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):  # a Windows console may not be able to show a snippet; never crash a report over it
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    limits = dict(DEFAULT_LIMITS)
    try:
        for item in a.limit:
            k, v = item.split("=", 1)
            limits[k] = int(v)
        problems, quotes = [], []
        for d in a.drafts:
            if not a.quotes_only:
                for hit in check_document(d, a.sources, a.max_run, limits):
                    problems.append((d, hit))
            for hit in find_quotations(Path(d).read_text(encoding="utf-8", errors="replace"), a.max_quote):
                quotes.append((d, hit))
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}")
        return 2
    for d, h in problems:
        print(f"{d}: {h['words']} words match {h['source']} (limit {h['limit']}), from word {h['doc_word_index']}: \"{h['snippet']} ...\"")
    for d, h in quotes:
        print(f"{d}: quotation of {h['words']} words (limit {a.max_quote}): \"{h['snippet']} ...\"")
    if problems or quotes:
        print(f"{len(problems)} run(s) over the copying limit and {len(quotes)} long quotation(s). Paraphrase, or get the publisher's permission.")
        return 1
    print(f"OK: no run over the limit in {len(a.drafts)} draft(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
