# Read and discard: grounding on sources the lab may cite but not copy

Status: proposed, written 2026-10-09. Nothing here is built or admitted. Untested design.
Provenance: drafted by `claude-opus-5-5` from a read of `sources/catalog.yaml`, `sources/README.md`, the cf-evidence-loop code and the lab's local-model scripts; reviewed by the session controller. No web page was fetched for it.

## The problem

Some sources may be linked but not copied (`redistribute: false`, or terms that bar copying). We want to check one claim against such a page without keeping the page. Reading a page into a model is still access, so this mode never widens what `access` allows.

## The rule

- Only `access: fetch` or `api` sources can ever be read this way, and only when a person has set `reference_use: read_discard` or `quote_short` for them. The default is `none`.
- `manual`, `request` and `forbidden` sources, and any host whose robots.txt disallows us (AI-agent groups included), is unreadable or unknown, stay link only. A person reads them in a browser.
- Where the terms are silent, the answer is "not permitted until a person decides". The first step for a source we want is to ask its publisher (see `2026-10-09-permission-and-review-letters.md`).
- A person pasting a page into a model is that person's choice under the publisher's terms. CFTR2's terms bar republishing "any portion": keep the link and our own words, never a quote, until written permission arrives.

## What "discard" means, honestly

Anything fetched with an agent's web tool or a shell command is stored in the agent session's transcript on disk (PDFs in full), so those routes are not read-and-discard. In the design below the page text never enters the agent session and is never written to disk by the reader: it lives in the reader's memory and the local model's memory, then is dropped. Kept: the URL, access date, SHA-256 of the bytes, the verdict, at most one checked short quote, and a discard record.

## The reader

A program (not a model):

1. Checks the catalog before any connection: the source exists, `access` is fetch or api, `reference_use` allows it, and the host matches.
2. Fetches within the existing network rules (pacing, robots.txt with AI-agent groups), with a byte cap; stops on 402, 403, 429 or an AI-use refusal (a `TDM-Reservation: 1` header or a `noai` tag).
3. Sends the text inside an `<untrusted_page>` boundary straight from memory to a local model with no tools and a JSON schema: verdict (supports, contradicts, not found in the text read), at most one quote, a location hint. No prompt files, no thinking files.
4. Checks the reply in code before dropping the page: the quote is an exact substring, within the word cap (8 by default, 0 under `read_discard`), a complete sentence, and holds the claim's key term or number. If any check fails, no quote is kept.
5. Prints one JSON record, then drops the page. The verdict is labelled inference; only a checked quote counts as T1.

It can ground one named, checkable claim on one page. It cannot ground absence beyond "not found in the N characters read", tables or figures, claims spanning pages, or anything from a source that is not fetch or api.

## Catalog change

New fields: `reference_use` (`none` | `link_only` | `read_discard` | `quote_short`, default `none`), `reference_decided` (date a person decided), `reference_basis` (the sentence of the terms that allows it), optional `max_quote_words` (default 8). A checker refuses `read_discard` or `quote_short` unless access is fetch or api, robots allows, `terms_url` is set, `claim_label` is verified, and both `reference_decided` and `reference_basis` are set. The evidence tool gets a matching gate beside its existing one.

## Tests before use

1. Sentinel: a local fixture page with a unique string; afterwards the string is in no file, log or output the lab's tools or Ollama write (operating-system level traces are out of this test's scope).
2. A made-up quote is refused; a 9-word quote is refused under the default cap.
3. An injection fixture cannot change the schema or get a quote through that is not on the page.
4. The gate refuses every manual, request and forbidden entry, and any entry left at `none`.

## Related findings (2026-10-09 log review)

- An agent's web tool does not check the catalog: on 2026-10-07 two requests went to a `manual` and a `forbidden` host; neither returned content. A hook that refuses such hosts is proposed to the maintainer (a settings change, human only), and should be in place before the reader is built.
- `delegate.ps1` leaves prompt and reply files in a temporary folder; the evidence tool keeps up to 700 characters of openFDA label text in a record. Both are noted for their owners.

## Decisions asked of the maintainer

1. Accept the rule that read-and-discard never applies to manual, request or forbidden sources?
2. Which fetch or api sources, if any, get `read_discard` or `quote_short`, after reading each source's terms?
3. Ask publishers for permission first (CFTR2 letter is drafted, unsent)?
4. Should the openFDA label excerpt shrink to the same quote cap?
5. Add a hook that blocks the agent web tool for hosts the catalog marks manual, request or forbidden?
6. Which repo owns the reader: cf-skills (beside the evidence tool) or cf-research?
