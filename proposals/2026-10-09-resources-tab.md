# A Resources tab: where people with CF, families and care teams can find help, by country

Status: proposed, written 2026-10-09. Nothing here is built or admitted. Untested design.
Provenance: drafted by `claude-sonnet-5-5` from a read of `sources/catalog.yaml`, `sources/README.md` and `proposals/2026-10-09-read-and-discard.md`; revised after a blind review by three `claude-opus-5-5` reviewers (compliance, community, engineering; findings were leads, checked by the controller against the catalog and the lab's records; no local model took part). No web page was fetched. The lab has not read any organisation's resource page, so this file makes no claim about what any organisation offers.

## What is asked

A tab on the site that gathers, in one place, the help people can use, by country, with a world view and a global section for things anyone can use. First the Cystic Fibrosis Foundation's resources, then four more countries. Grouped by who it is for (people with CF, families and carers, care teams) and by what it helps with.

## What the lab knows today (observed in the catalog, 2026-10-09)

- `www.cff.org` is `access: manual` in all twelve Foundation entries (five annual reports, the technical supplement, a COVID handout, five highlights); robots.txt was unreadable on 2026-10-05 and the notes say the site "refuses automated clients". An agent never fetches it, and the read-and-discard design (proposed, not built) excludes manual hosts. A person reads those pages in a browser.
- Other catalogued bodies: UK (Cystic Fibrosis Trust, `fetch`), Australia (Cystic Fibrosis Australia, registry run by Monash University, `fetch`), Ireland (the Cystic Fibrosis Registry of Ireland at University College Dublin, `fetch`; this is a registry, not a support organisation), Canada (Cystic Fibrosis Canada, `manual`), Europe-wide (ECFS Patient Registry on `pr.ecfs.eu`, `fetch`; the society's main site is not in the catalog). France's entry is a variant list, not a registry report.
- Each of those entries covers one document (a report or a list). The terms on record are for those documents. No organisation's resource pages have a catalog entry, a robots check or a terms reading yet.
- Not known: whether any Foundation program serves people outside the United States. The lab will not state that from memory; each row records what its own page says.
- The drafter's inference, not in the catalog: that all five proposed countries publish in English. The lab's rule that it does not machine-translate health information is policy, not a catalog fact.

## Design

### 1. Data first

A new folder `resources/` in cf-research holds `resources.yaml` (starting with `schema_version: 1`) and the checker. One row per resource:

| Field | Meaning |
| --- | --- |
| `id`, `name`, `organisation` | Stable slug (also the page anchor); the resource's name as the organisation writes it; who runs it (spelled from one allowed list) |
| `countries` | A list of ISO 3166-1 alpha-2 codes (`GB`, not `UK`), or the region values `europe` and `global` |
| `url` | One https link to the organisation's own page |
| `audience` | One or more of `adult_with_cf`, `teen_young_adult`, `parent_carer_of_child`, `partner_family_carer_of_adult`, `care_team`, `everyone` |
| `category` | One or more of the headings below |
| `summary` | One sentence in the lab's own words, at most 25 words, no quote marks, no digits other than a year, no URL |
| `who_for` | A fixed value: `stated_anyone`, `stated_countries`, `stated_region`, `not_stated`. Never the organisation's wording. The page shows a fixed line, not the row's text |
| `format` | `online`, `in_person` or `both` |
| `lang` | BCP 47 language of the linked page |
| `source_id` | A catalog entry of the new kind `resource_pages` whose landing-page host equals the host of `url`. An entry for a registry report or a highlights file fails this rule |
| `basis` | `person_read_in_browser`, or `script_fetch` (only for a host whose entry has `access: fetch`, a robots check within six months and the terms read) |
| `summary_attested_by`, `summary_attested_on` | A person states the summary is their own words and matches the page. No script can check this: the page text is never kept, and manual hosts are never fetched. It is an attestation, not a test |
| `confirmed_by`, `confirmed_on`, `content_hash` | A person's confirmation of the row as it stands. The hash covers name, url, summary and who_for; a change after `confirmed_on` makes the row unconfirmed. A later script refresh cannot overwrite it |
| `last_checked`, `status` | When a person last opened the page; `active` or `retired` (with a date) |

The file holds no logos, no organisation text beyond the resource's name, no phone numbers or emails (links only; numbers go out of date and the lab writes none from memory), no personal stories, and no outbound-click logging.

The checker (`tools/resources/check_resources.py`, standard library only, with a `--today` flag so tests are deterministic, run in check-all and in CI) enforces what a script can see: schema and enumerations, https, the `source_id` rule, the word cap, no phone or email, valid ISO codes, hashes. It prints a "due for a person to re-check" list.

Freshness. A row is due for a re-check at 180 days (90 for the categories Money, insurance and benefits and Medicines and access, where details change); the checker warns and lists it. At 365 days (180 for those two categories) the page marks the row "not re-checked since <date>". CI never fails on a date. No tool, test or link checker requests any `url` whose catalog `access` is manual, request or forbidden; those links are re-checked by a person from the due list.

### 2. Who it is for, and what it helps with

Audiences: adults with CF; teens and young adults with CF; parents and carers of children with CF; partners, family and carers of adults with CF; care teams; everyone.

Categories, in this order (the page shows only those with confirmed rows):

1. Urgent and crisis support
2. Newly diagnosed and learning the basics
3. Find a care team or clinic
4. Mental health and emotional support
5. Medicines and access
6. Money, insurance and benefits
7. Daily life: school, work and travel
8. Growing up and moving to adult care
9. Relationships, fertility and family planning
10. Transplant and advanced lung disease
11. Grief and bereavement
12. Community and peer support
13. Clinical trials, research and registries
14. Advocacy and getting involved

Only public pages run by a named organisation are listed, not private groups, social-media pages or forums. Any infection-control wording about meeting other people with CF comes from a person-checked source, not memory.

### 3. Which countries first

Five: Australia, Canada, Ireland, the United Kingdom and the United States (alphabetical on the page). Reasons: each has a registry report already in the catalog, so the lab knows the field there; the drafter infers that all publish in English, so the lab can check every row's wording itself. The organisations that run resource pages, and their terms, still need their own entries (for Ireland the catalogued body is the registry, so a support organisation must be found and entered). A `europe` and a `global` section hold bodies that span countries, once each has an entry for its main site. France is a natural sixth when a French reader is found. The maintainer may swap any country.

The page says plainly: five countries are listed so far; the others are not listed yet, which does not mean nothing exists there; your CF care team or your national CF organisation is a good first stop. It gives no prevalence or count of countries.

### 4. How the rows are made

| Host `access` | How |
| --- | --- |
| `manual` (the Foundation, CF Canada) | The maintainer opens the pages in a browser and sends only names, links and a one-line note in their own words. No agent reads a saved copy of the page, fetches the host, or uses a headless browser. An agent may structure what the person sent; a person writes or approves every summary against the page |
| `fetch` (after a `resource_pages` entry exists, with robots and terms read by a person) | A script requests the single page named by each row, one at a time, never crawling a directory; keeps only the name and URL; drops the page text. A person writes the summary |

Before any row: the terms of use of the host are read by a person (robots.txt alone is not enough), and the organisation is written to (the lab's permission-letter route), so no organisation first learns of the listing from the page. Rows are confirmed by a person; this is a script plus a person, not tier T3, because no second model sees the page.

### 5. The tab

- A page `/resources/`, generated by a tool in the site repo from a pinned snapshot of `resources.yaml`: `data/resources.snapshot.yaml` and a lock file with the cf-research commit and the file's SHA-256. The generator reads only the snapshot and refuses a `schema_version` it does not know. A small tool moves the pin, in its own pull request. (Reading cf-research's main would repeat the merge-order failures already seen.)
- It opens with one line, then the country chooser (a plain alphabetical list with "Global" and "My country is not listed"). The Foundation's rows sit under United States like any other country's. The page has no group "first".
- Opening line: "This page gathers help that organisations in a few countries offer, so you can start from one place." Directly under it, with no number in it: "If you or someone with you is in danger or in crisis now, contact your local emergency number or your CF care team's urgent line. This page is not a crisis service." (This follows the lab's current rule; the open board row on the emergency-number wording may change it.)
- Each row shows the name, the summary, a fixed "who it is for" line ("Check the organisation's page for who can use this"; for `not_stated`: "The organisation's page does not say who can use this. Check with them before relying on it"), the format, and the date last checked.
- Note: "We list resources so they are easier to find. Being listed is not a recommendation, and the order is not a ranking. We are not part of any organisation here and are not paid by them. Each organisation decides who it can help, so its own page has the last word. Information changes: check with the organisation and your care team. This is research, not medical advice, and not advice about money or medicines for your situation. Spotted a broken link or a mistake? Tell us on the Contact page." Listing rule: a named organisation, a public page, free to view. Organisation names only, nominative: no logos, colours or wording that suggest partnership; no organisation name in the tab's title. (The "65 Roses" phrase is the Foundation's registered trademark per the Foundation's 65 Roses page as recorded in `cf-research-context`, checked 2026-10-03; the page does not use it.)
- Filters by audience and category are written into the static HTML and work with CSS only (`:has()` over checkboxes and server-written `data-` attributes), so there is no layout shift and no new script; with CSS unsupported every group is visible. No analytics or tracking.
- The generator computes the page's bytes and fails at 80% of the site's cap. When a page passes that, `/resources/` becomes an index with `/resources/<code>/` per country. A country with no confirmed rows shows a plain "none listed yet" line.
- Tests: page weight, layout shift, scripts-off shows everything, each filter's row count, Chromium and WebKit, headings (one level-2 per country), and that no outbound link host is manual, request or forbidden in any automated test.

### 6. The map (later, optional)

Only after the list ships. Boundary data from a public-domain source whose licence is read and recorded in the catalog, with its URL and SHA-256 (a third-party download needs the maintainer's say-so; a note on disputed borders). A script simplifies it into an SVG under a size budget (25 KB suggested), regenerated and diffed in CI. Shapes are decorative; each listed country is one link around a marker of at least 44 by 44 CSS pixels, with a callout for small countries such as Ireland, a visible focus ring, a title, and a text equivalent that is the country chooser. Uncovered countries are neutral ("not listed yet"), never shaded as "no help". The page works without it.

### 7. What this does not do

No ratings, rankings or "best of". No medical advice, eligibility advice or advice about money. No copying of an organisation's text or logo. No private groups. No claim about how well any resource works. No row from memory; a model may structure what a person sent, but a person writes or approves every summary. No fetch of any manual, request or forbidden host by any tool.

## Order of work (the smallest slice first)

1. Prerequisites: the web-tool block for catalog hosts marked manual, request or forbidden (board T-0108), so no agent request can reach cff.org.
2. cf-research: the schema, the checker and its fixture tests (failing and passing rows), wired into check-all and CI; a `resource_pages` kind in the catalog and one entry for the Foundation's resource pages (`manual`).
3. The maintainer reads about five Foundation pages and sends names, links and notes; an agent structures the rows; the checker runs; the maintainer attests and confirms. This also measures the minutes per row for a re-check.
4. Site repo: the snapshot, the lock, the generator, the list-only page, and its tests.
5. The other four countries, one at a time, each with its catalog entry, terms reading and letter first.
6. The filter, then `europe` and `global`, then the map.

## Decisions asked of the maintainer

1. The five countries (Australia, Canada, Ireland, UK, US), or your own five?
2. Will you read the Foundation's pages in your browser and send names, links and one-line notes? (About five pages for the first slice; the cost per row is measured then.)
3. The tab's name and place in the menu ("Resources", after "CF"?).
4. Write to each organisation before listing it (recommended)? Who signs, given the open identity decision on board T-0117?
5. The freshness numbers (180 and 365 days; 90 and 180 for money and medicines): keep or change? Who re-checks each country?
6. The urgent-help line: keep the "local emergency number" wording until the emergency-number rule (board T-0104) is decided?
7. A self-hosted public-domain world map, later: acceptable as a third-party download?
