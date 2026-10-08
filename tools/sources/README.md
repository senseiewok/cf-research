# tools/sources

Fetches the source documents the repo is permitted to fetch, and prints instructions for the ones a human has to download.

Reads [`sources/catalog.yaml`](../../sources/catalog.yaml). Writes `sources/manifest.json`. Does not redistribute anything.

## Install

```bash
pip install -r requirements.txt
```

Python 3.9+.

## Run

```bash
python fetch_sources.py              # fetch permitted sources, then list what you must download
python fetch_sources.py --manual     # print the browser download list only
python fetch_sources.py --verify     # re-hash local files against the manifest
python fetch_sources.py --only ecfspr-adr-2022 --verify   # check just the named entries
python fetch_sources.py --record     # checksum files already on disk, e.g. ones you saved by hand (no network)
python fetch_sources.py --only ecfspr-adr-2022
python fetch_sources.py --dest /path/outside/the/repo
python fetch_sources.py --force      # re-download files that already exist
```

Default download directory is `sources/downloads/`, which is gitignored.

## What it will not do

It downloads an entry only when the catalog says `access: fetch` or `access: api`. Entries marked `manual`, `request`, or `forbidden` are reported, never retrieved — and there is deliberately no user-agent option, no browser-automation path, and no retry-with-different-client behaviour. If a publisher refuses automated access, the answer is to ask them or to download by hand.

Every download goes through the `Downloader`, which follows the evidence skill's network rules (`cf-skills/.claude/skills/cf-evidence-loop/references/NETWORK-RULES.md`) and imports its robots rule rather than keeping a copy:

- it identifies itself with the project user agent and no other, over https only;
- it reads the host's `robots.txt` once and obeys it, including groups that name AI agents (this tool is one); an unreadable `robots.txt` means no download;
- it follows redirects only within the same host and checks `robots.txt` again for the new path;
- it makes one request per file and never retries; it pauses at least two seconds between requests to a host, longer for a `Crawl-delay` (up to 30 seconds); it refuses a file over 200 MB and leaves no partial file.

The evidence skill must be next to this repository (`../cf-skills`) or named by `EVIDENCE_SKILL_DIR`, and `requests` installed; without them the tool refuses to download and `--manual` still works. `python tools/sources/fetch_sources.py --audit-robots` reads each fetch host's `robots.txt` once, requests no document, and reports whether every `fetch` entry is still allowed under the current rule (exit 1 if any is refused); run it before trusting a `robots: allow` field in the catalog. The conduct tests are `python tools/sources/test_fetch_conduct.py`.

`--verify` is the part that matters for reproducibility. The manifest records the SHA-256 of the exact bytes an extraction ran against, so a stranger can download the same report from the publisher, confirm the hash, and reproduce the result without us hosting a single copyrighted page.

## Exit codes

`0` success. `1` a download failed, or `--verify` found a changed or missing file, or `--verify --only` named an entry that has no recorded hash yet, or `--verify` could check no file at all (every entry `NO HASH`, or an empty manifest), or the existing manifest could not be read (it is then left as it was). A full `--verify` lists every entry with no recorded hash as `NO HASH` (save the file and run `--record`) but still exits `0` when the files it could check all match. `2` unknown catalog id, including `--verify --only` with an id the manifest does not have.

## Tests

```bash
python -m unittest -v     # every test in this folder; no network, no real catalog
```

Run it from `tools/sources` with PyYAML installed (`pip install -r requirements.txt`); the tests are not standard library only. Each `test_*.py` also runs on its own, for example `python test_check_numbers.py`, and the four checker tests (`test_check_catalog_fields.py`, `test_check_numbers.py`, `test_check_source_overlap.py`, `test_make_manual_downloads.py`) then accept an optional path to a candidate script to test instead of the real one, and print `VERIFIED`.

`check_catalog_fields.py` checks that every key used in `sources/catalog.yaml` is documented in the Fields table of `sources/README.md` (exit 1 and the key names when one is not; a documented key nothing uses is only a warning). `python test_check_catalog_fields.py` tests it on fixtures.

## Guards against invented text

Three small checks, each with a test that fails when the check is broken. None of them decides whether a claim is true; they stop the failures a model makes while writing a note from reports.

```bash
python tools/sources/check_source_overlap.py DRAFT.md            # copying: runs of words shared with sources/downloads; also quotations longer than 8 words
python tools/sources/check_source_overlap.py --quotes-only DRAFT.md   # the quotation check alone; needs no downloads
python tools/sources/check_numbers.py --evidence EVIDENCE.txt --allow-file ALLOW.txt DRAFT.md
```

The lab's rule is to paraphrase and not to quote a source's sentences in a public note; names, titles, labels and field names are short enough to pass. `check_numbers.py` flags every number in a draft that the evidence does not contain. Use narrow evidence, the quotations and summaries the note was written from, because a whole folder of reports contains almost every number and the check then proves nothing. A number that is arithmetic goes in the allow file with its sum (`564  # 127 + 437`), and a line with no reason is an error. A pass means the numbers exist in the evidence, not that each sits on the right claim; a second model and a person still read the claims. The tests are `test_check_source_overlap.py` and `test_check_numbers.py`.

## Saving a manual document

A source marked `manual` is downloaded by a person in a browser. Save it into `sources/downloads/` under the filename `--manual` shows, run `--record` to checksum it, then `--verify`. The manifest never stores a machine path. A `--record` or download run keeps the hash an earlier run recorded for a file that is not on this machine (marked `present: false`), so recording on a computer that holds only some of the files does not erase the others' hashes; it never makes one up.

## The manual downloads page

`sources/manual-downloads.md` lists every source a person must download by hand (`manual`), the ones to ask for (`request`) and the ones not to download (`forbidden`), with links, the file name to save each one under, and a prompt to paste into an agent when the files are saved. It is generated from the catalog and must not be edited by hand:

```bash
python tools/sources/make_manual_downloads.py           # write the page
python tools/sources/make_manual_downloads.py --check   # exit 1 if it is out of date
python tools/sources/test_make_manual_downloads.py      # prints VERIFIED
```

Run it after any catalog change that touches a `manual`, `request` or `forbidden` entry.

## Checking a draft against the sources

`check_source_overlap.py` measures how many consecutive words a draft shares with the documents in `sources/downloads/` (PDFs through `pdftotext`, or `.txt` files). Most reports may be read, cited and linked but not reproduced, and the UK report asks for permission before its content is reproduced in publications, so a public note should paraphrase and quote only briefly:

```bash
python tools/sources/check_source_overlap.py DRAFT.md [DRAFT.md ...]
python tools/sources/check_source_overlap.py --max-run 20 --limit CFFPR=15 DRAFT.md
python tools/sources/test_check_source_overlap.py          # prints VERIFIED
```

The default limit is 25 words in a row, and 10 for files whose name starts with `UKCFR`. Case, punctuation and line breaks are ignored, and the report prints the first six words of a run, not the run. Exit 0: nothing over the limit; 1: a run over it; 2: bad input. It finds copying; whether a shorter quotation is fair is a person's call, and it does not replace the catalog's `redistribute` and terms fields.
