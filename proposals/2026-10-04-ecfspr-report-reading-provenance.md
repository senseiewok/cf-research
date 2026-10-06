# Provenance: how the ecfspr-report-reading skill was made

Status: **placed in `.claude/skills/ecfspr-report-reading/` of this repo at the human's direction (2026-10-04); content review is part of the PR review.** It is not in `hq/.claude/skills/`. A skill under `.claude/` is agent instruction, so a human reads the diff before it is relied on.

### Inputs

- Five ECFSPR annual reports (2020 to 2024), fetched by `tools/sources/fetch_sources.py` from `pr.ecfs.eu`, which the catalog permits. They are git-ignored and checksummed in `sources/manifest.json`. Text came from `pdftotext -layout`.
- Nothing from `www.cff.org`, `databases.lovd.nl` or `cftr2.org`. No browser automation was used.

### Steps

1. A script (not a model) chose 20 pages: contents, introduction and general considerations, and the summary table, for each report. It tolerates the 2022 report's corrupted letters in headings.
2. Local `qwen3.8:27b-64k`, fast profile, no tools, input inside `<untrusted_page>` tags, two samples, was asked for facts each backed by an exact quote. A script kept a fact only if its quote appears verbatim on that page. First run: 17 and 40 facts returned, 47 unique facts kept after checking.
3. A script built a contents table and a printed-versus-PDF page offset for each year.
4. The same local model drafted skill text from those verified items in two samples.
5. **Review found both model drafts wrong in the same place.** Each said the 2024 definition of chronic *Pseudomonas aeruginosa* infection "returns to mentioning only the modified Leeds criteria", and that 2020 mentioned only those criteria. The source sentences in all five reports include the antibody criterion too. Cause: the quote supplied for one year was a truncated prefix of its sentence, and the model built an "only" on it. The mechanical checks (citation keys exist, numbers appear in the material, quotes are verbatim) all passed. They cannot detect a false inference.
6. The final text in `SKILL.md` was therefore rewritten by the reviewing agent (`claude-sonnet-5-5`), keeping the model drafts' section structure and using only statements checked against the full sentences. A script then confirmed each quoted phrase occurs verbatim in the named year's report and each layout statement matches the contents tables (the final counts are in step 7). The real cross-year differences it records (where the definitions list lives, a stale page cross-reference in the 2021 report, the 2024 three-year footnote) came from that review, not from the model.

7. **Blind review by a frontier model.** One read-only subagent (the Agent tool's `fable` model override, a different model from the author, same provider) was given the skill and the five extracted texts, told to use no outside knowledge, and asked to check every statement. It had no tools beyond reading and no access to the earlier findings. It found three overbroad statements. A script confirmed all three against the source before any edit:
   - "The 2020 contents have none of these" (late sections) was wrong in effect. The 2020 and 2021 contents list an unnumbered "Data Quality" entry. My contents parser captured only numbered entries, so my own check had a blind spot.
   - The population section described two footnotes. Every report has three, and the United Kingdom note belongs to the third.
   - The "three years" quote was placed under the United Kingdom bullet. It is in the first footnote.
   The text was corrected, and the verifier now also checks the raw contents pages and the footnote structure of every year (35 quotations, 30 layout statements and 18 raw-contents and footnote checks, all passing). The reviewer found no quotation false.

### What this does not cover

- Only the contents, introduction and summary table pages were read for each year. The data sections, appendices and the CFF Technical Supplement were not.
- The printed-versus-PDF page offset was checked at one heading per report.
- The skill was not tested in any agent client, and its description has not been evaluated for triggering.
- This is one small run on one registry. It does not show that the local model drafts skills reliably, or that a frontier reviewer will always catch what scripts miss. The review used one model once.

### Human review items

- Read `.claude/skills/ecfspr-report-reading/SKILL.md` in full. Placement in `research` was decided by the human; content has not been approved line by line.
- Decide whether the `cf-research-context` "Registry numbers" section and this skill overlap enough to merge.
- The 2022 report's corrupted letters are described but not normalised. A deterministic repair rule has not been written.
