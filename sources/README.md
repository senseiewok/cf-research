# Sources

[`catalog.yaml`](catalog.yaml) enumerates every external source this repo depends on and records, per source, **what we are permitted to do with it**. Tools and agents read the catalog instead of hardcoding URLs or assuming terms.

Nothing here is a legal opinion. The catalog records what we tested and when; where we have not tested, it says so.

## Why a catalog rather than links in prose

Three things kept going wrong without one. URLs drifted and nobody noticed. The same source got re-checked by different people who did not know it had already been checked. And the question "are we allowed to download this?" was answered from memory, which is exactly the kind of thing that should not be answered from memory.

## The one rule

**`access` is a permission, not a hint. An agent may automate retrieval if and only if `access` is `fetch` or `api`.**

| `access` | Meaning | What an agent may do |
| --- | --- | --- |
| `fetch` | Automated download is permitted | Download the file |
| `api` | Queryable API, permitted | Call the API, honouring the API's published rate limits and any `etiquette` |
| `manual` | A human must download it in a browser | Print instructions. Never fetch, and never use a headless browser to make the request look human |
| `request` | Needs an application or written permission | Stop and tell the human what to ask for |
| `forbidden` | The publisher has explicitly refused automated access | Never touch it. Link to it in prose only |

If a source is not in the catalog, it has no permission. Add an entry first, with evidence.

## Fields

| Field | Meaning |
| --- | --- |
| `id` | Stable slug. Referenced by extracted rows, so never renamed once used |
| `title`, `publisher`, `kind` | What it is and who publishes it |
| `data_year` | For annual reports, the year the data describe, **taken from the title page** |
| `access` | The permission above |
| `robots` | `allow`, `disallow`, `unknown`, or `unreadable` (we could not retrieve robots.txt) |
| `robots_checked` | Date we last looked. Treat anything older than six months as stale |
| `url` | Direct file URL, or `null` when only the landing page is known |
| `landing_page` | Where a human goes to find it |
| `filename` | The local name the fetcher writes |
| `redistribute` | May we commit or host the file itself? Set once under `defaults:` (`false`), not per entry |
| `citation` | The publisher's own preferred citation string, where they state one |
| `claim_label` | `verified` / `unverified` / `hypothesis`, same convention as the rest of the repo |
| `warning` | A machine-readable gotcha, e.g. `filename_year_mismatch` |
| `notes` | Anything a future reader needs, including what was *not* checked |
| `pages` | Page count of the PDF, when someone has checked it. Recorded for one report so far |
| `text_layer` | `embedded` when the PDF has selectable text. No entry needs OCR yet |
| `priority` | Optional hint for which documents to obtain first. Only `high` is used so far |
| `reference_doi` | DOI of the paper that describes the source, where one exists |
| `base_url` | Root URL of an API, for entries with `access: api` |
| `api_root` | Informational API root for an `api` source that is under admission review. The evidence client reads only `base_url`, so a source with `api_root` and no `base_url` cannot be reached by it |
| `terms_url` | The page that states the API's usage policy. Required for every `api` source: the evidence client refuses to call one without it |
| `max_rps` | Our self-imposed ceiling in requests per second for that source. The client never exceeds 3 per second per host whatever this says |
| `published_limit` | What the publisher's page says about limits, or plainly that none was stated and the date it was read |
| `auth` | `none`, `optional_api_key` or `api_key_required`: whether the API needs a key |
| `etiquette` | Free text: any politeness rule the API asks of clients |

## What is and is not stored here

**In the repo:** the catalog, the generated manifest, and extracted values with provenance.

**Not in the repo:** the source PDFs. They are third-party copyrighted publications. We link to them, record a SHA-256 of the exact bytes we used, and leave distribution to the publishers. Downloads land in `sources/downloads/`, which is ignored by git.

A checksum is what makes an extraction reproducible without redistributing anything: a stranger downloads the same report, confirms the hash matches, and re-runs the pipeline.

## Usage

```bash
pip install -r ../tools/sources/requirements.txt
python ../tools/sources/fetch_sources.py            # fetch what is permitted, list the rest
python ../tools/sources/fetch_sources.py --verify   # re-hash local files against the manifest
python ../tools/sources/fetch_sources.py --manual   # print browser download instructions only
```

See [`../tools/sources/README.md`](../tools/sources/README.md) for options, and [`../landscape/registry-data-access.md`](../landscape/registry-data-access.md) for the access testing this catalog is built on.

## Adding a source

1. Check `robots.txt` and the publisher's terms. Record the date.
2. Add an entry. Set `access` from what you found, not from what would be convenient.
3. If you could not check something, set `claim_label: unverified` and say so in `notes`.
4. If a publisher refuses automated access, set `access: forbidden` and leave it refused. Do not look for a way around it.
