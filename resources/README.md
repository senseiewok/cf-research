# Resources

The list behind the site's Resources tab: help that organisations offer to people with CF, families and care teams, one row each, **links only**. The design and the reasons for each rule are in [`../proposals/2026-10-09-resources-tab.md`](../proposals/2026-10-09-resources-tab.md).

| File | What it is |
| --- | --- |
| [`resources.yaml`](resources.yaml) | One row per resource. The only data file. |
| `../tools/resources/check_resources.py` | Checks the rows against the schema and the source catalog (`python tools/resources/check_resources.py`). |

## What a row is

A name, one link, who it is for, what it helps with, a one-sentence summary in the lab's own words, and the dates a person last looked. It holds no phone number, email address, logo, copied text, price, count, staff name, individual study or campaign.

## How a row is made

1. The organisation's host has an entry of kind `resource_pages` in [`../sources/catalog.yaml`](../sources/catalog.yaml) that says what the lab may do with it. For a `manual` host (the Foundation's) no tool or agent requests the page: a **person opens it in a browser**.
2. The person sends the name, the link, who the page says it is for (or "not stated"), the date, and a note. An agent may structure what the person sent. The summary is the lab's own words; the person attests it matches the page.
3. `python tools/resources/check_resources.py --print-hashes` prints each row's hash; the person's confirmation (`confirmed_by`, `confirmed_on`, `content_hash`) goes on the row. Changing the name, link, summary or `who_for` afterwards makes the hash wrong, so the row counts as unconfirmed until a person confirms it again.
4. A name seen on a slide or in an overview, and not yet opened by a person, is a `lead_only` row at most: it is never confirmed or published.

## What the checker does and does not do

It checks what a script can see: fields, allowed values, https, the catalog link, the word cap, no phone or email, valid dates, hashes, and which rows are due for a person to re-check (180 days; 90 for money and medicines). It **cannot** check that a summary is not copied or that a page says what the row says: that is the person's attestation. It never requests any link.
