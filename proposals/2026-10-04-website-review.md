# Review: the lab website against the lab's own rules

Status: review note. Read 2026-10-04 from the website's source (not the deployed page). The site's public address is senseiewok.ai; this note does not record how or where it is built or deployed (that stays in the private site repo). Each finding is a task on `hq/tasks/BOARD.md`; ids given. Standard applied: `hq/README.md`, `cf-research-context`, `lab-voice`, `security-browsing`.

## Verdict first

The site is in much better shape than the YouTube About text. It carries "research, not medical advice" four times, says "independent" three times, dates its checks, labels the population figures as estimates, declines to list dosing or eligibility, and separates the lab's scope from the Foundation's mission. The modulator history paragraph is a model of how to write about drugs without advising. Most of what follows is tightening, not repair.

## Findings

| # | Finding | Rule | Task |
| --- | --- | --- | --- |
| 1 | `hq/LICENSE.md` names no license, has a placeholder contact, and `research`, the website repo, and `skills` have none. The site says "open source" eleven times | Without a named license, default copyright applies; "open source" is then untrue | T-0001 (P0) |
| 2 | "M365 Copilot drives the lab's task scheduling, turning intentions into a calendar of bounded, reviewable next steps" | There was no task list in `hq` until today. The sentence describes something that did not exist | T-0015 |
| 3 | The 65-rose field is JS-rendered (`<div id="rose-field">` is empty in source). README and `cf-research-context` promise exactly 65 and say to count rendered instances | Only a browser can verify; nothing currently does | T-0016 |
| 4 | Google Fonts and an `images.unsplash.com` preconnect. No Unsplash image is used anywhere | Every visitor's IP goes to two third parties. For a CF-community site this is an avoidable leak; `security-browsing` section 2 spirit | T-0017 |
| 5 | "Independent" appears three times; "not affiliated" never, while the donate call and CFF history dominate the page | `cf-research-context`: "Do not ... imply affiliation." One explicit line closes it | T-0018 |
| 6 | Tagline "CF - until there's a cure" | `lab-voice`: "Hope without a forecast." Not a claim the lab cures anything; a human call whether to keep it deliberately | T-0014 |
| 7 | Field guide names three locally installed models; `hq` has profile skills for only some of them | Public claims about the cast should match the skills that govern it | T-0019 |
| 8 | Approval years (1993, 2019, 2024) cite a CFF history page and one DailyMed label | Accurate as far as I know, but a per-drug FDA or DailyMed citation is the standard the page sets for itself elsewhere | T-0023 |
| 9 | A cosmetic mismatch in the site repo's internal contributor notes (details kept in that private repo) | Cosmetic; a contributor will trip on it | folded into T-0015 |

## What was checked and found fine

- `prefers-reduced-motion` is handled in CSS.
- 51 ARIA attributes; skip link present; the rose field is `aria-hidden`.
- The site's deployment practice was read against `security-baseline` and nothing was found; the details are deliberately not recorded in this public note.
- The population figures (40,000 US, 105,000 worldwide) are labelled as a 2022 estimate with a check date, not presented as current counts.
- No patient data, no dosing, no eligibility statements, no logos.

## Not checked

- The deployed page itself. Only the source was read; a deploy could differ.
- Lighthouse / accessibility beyond ARIA presence.
- A voice skill in the website repo was not compared against `lab-voice` in `hq`. Two voice skills in two repos is a likely drift point worth a task if they are both maintained.
