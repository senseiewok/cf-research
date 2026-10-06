# Security review of the 2026-10-05 work: misses found, fixes made, and what is left

Status: review, written 2026-10-05 from the week's commits and the way the sources were read. Scope: the network and terminal conduct of the evidence loop, the way third-party pages were read, and the review process itself. It is a self-review by the agent that did the work, plus one local-worker review of the T-0052 diff; the full tier's blind challenger from another model family was not available, so this does not clear that row.

## Fixed in this batch (each has a test that failed before the fix)

| Miss | Fix | Where |
| --- | --- | --- |
| A response body of any size was read | A 5 MB cap, from a declared length or the streamed bytes (compressed bodies are counted after decompression) | cf-skills 129d43e (R10) |
| A JSON list, a missing result key or a non-numeric total became a crash or a silent "nothing found" | Error records that say what was wrong; a missing total is never read as zero | 129d43e |
| An error message could carry a whole response value, and a total of more than 4,300 digits crashed `int()` | Counts are at most 15 digits; the value shown is cut to 60 characters | 80dd0b1 |
| Publisher text could send terminal escape sequences to the reader | Escape sequences and control characters are removed before printing (timed on 5 MB adversarial inputs: 0.02 to 0.26 s) | 129d43e (R11) |
| `--contact` went into a header unchecked, and options given before the subcommand were silently dropped | A plain-email check; the options are no longer overwritten | 129d43e |
| robots.txt group matching used prefixes, let an empty agent match everyone and used only the first matching group | RFC 9309 matching; 320,000 random cases agree with an independent reference | 129d43e |
| A dry run did not list the requests it would make | Every planned request, robots.txt included, is listed | 129d43e (R9) |
| The page reader followed redirects to other hosts, used the system trust store, used the standard library's robots parser and ignored Crawl-delay | Same-site redirects only with robots re-checked, certificate checks through `requests`, the lab's own matcher, Crawl-delay honoured. Scratch copy only; see the open items | untracked reader script |

## Open, with a recommendation

| Item | Why it matters | Recommendation | Board |
| --- | --- | --- | --- |
| The reader is an untracked script | Reading terms pages is conduct code with no tests; two of the misses above were in it | A tested command with the evidence client's rules, for human-started reading of terms pages only | T-0066 |
| robots.txt groups that name AI agents do not apply to us | Cochrane's file names `Claude` as disallowed; the evidence package matches only its own token, so a catalogued fetch source on such a host would be fetched | A policy decision for a human: treat a disallow for the common AI-agent tokens as a refusal | T-0067 |
| The full tier's blind challenger and a human diff read are unticked for this week's commits | The gate says "unavailable means blocked"; the PRs say so | Say it in each PR (done) and let the human read the diffs; do not describe the tier as complete | none |
| The local worker's review row had never been run | The row exists to find what the controller missed; the one run found nothing true and the controller's own read found one gap | A one-command helper makes the row cheap; keep treating its output as a prompt | T-0068 |
| Terms and rate limits in the catalog were written by a model from pages it read | A wrong `claim_label: verified` would be a public false statement | A human reads the new entries; the PRs list them. The three `terms_url` values from T-0052 are still unconfirmed | T-0052 |
| An untracked test file from a stood-down session sits in the cf-skills working tree | An accidental `git add -A` would commit unreviewed code | Stage explicit paths (done); ask the owner to delete or hand it over | none |
| VS Code rewrote the tracked workspace file | An editor edit could ship by accident | Left uncommitted; the owner decides whether the setting belongs in the shipped file | none |
| `ukri-gtr` is now a source the skill may reach, with no provider | Harmless today (the client still enforces pacing, budget and the body cap); a provider needs its own shape checks | Add the provider with fixtures and checks | T-0069 |

## Things checked and found fine

- No change added a dependency, and the only new network code in a repo is the evidence client's own (already rule-bound). The reader is the one new network-facing script and it is not in a repo.
- No page text was given to the local worker or to any tool; the data boundary in `security-browsing` was not needed, and no instruction was found in a page. If that changes, wrap the text first.
- Certificate checks were never switched off. A self-signed certificate in the chain on one host (seen with the system trust store, not with the evidence environment's) was reported and worked around by using the evidence environment, not by disabling verification. It may mean this machine sits behind an inspecting proxy; worth knowing.
- Pages were read with robots.txt first, the site's Crawl-delay, one request each and an honest user agent. Where a site refused (www.cff.org, Cochrane, ClinGen) nothing was fetched and the catalog says why.
- The new scripts (`asciicanvas.py`, `check-ascii.py`, `delegate-batch.ps1`, `make_manual_downloads.py`) read files and write files only and open no network socket themselves; `delegate-batch.ps1` runs `delegate.ps1`, which talks to the local model on loopback. `delegate-batch.ps1` accepts item names with letters, digits, dot, underscore and hyphen only, so a name cannot walk out of its folder.
- A privacy scan ran clean on every commit.

## What this review cannot show

It is one agent checking its own work, with one weak second opinion. It does not show the absence of other misses. A blind review from a different model family, and a human read of the diffs, are still the way to find what neither of us looked for.
