# Review record: the CFTR model package (tools/cftr-model)

Status: record, written 2026-10-05 by the controlling agent from this session's work. It describes how the package was made and checked and what is still open. It is not a clearance: nobody with CF or a carer has read the page yet (task T-0081), and the page is labelled "Draft, in review".

## What was made

A browser model of the human CFTR protein from four published structures (6MSM, 5UAK, 8EIQ, 8EJ1), with research tools, and the generator, data builder and tests that keep every number on the page tied to the source files. The facts it rests on are in `landscape/cftr-structures.md`; the sources are catalogued in `sources/catalog.yaml` (five new entries, all `redistribute: true`).

## Who did what

| Role | Model | What it did | What it did not do |
| --- | --- | --- | --- |
| Controlling agent | Claude (this session) | Scoped the work, wrote the data builder, tests and page text, ran every check, owns the diff | Did not decide scientific wording alone: each changed claim was checked against a source line |
| Design, critique and renderer author | Fable, through the Agent tool | Storyboard, design critique, the WebGL2 renderer (written in an isolated copy and merged after tests), the navigation plan, a security and an SEO review of the site | Its findings were leads; each was confirmed on source before use |
| Blind reviewers | Opus (cloud), Qwen3.8 27B (local) | Read the captions, the page text and the plans without seeing each other's findings; flags were confirmed on source, never counted as votes | Neither approved anything for publication |
| Local worker as author | Qwen3.8 27B | Attempted to write the WebGL2 viewer: 0 of 6 attempts accepted by the tests (recorded in the model's profile skill in cf-lab) | Not used for the renderer |

## Checks that decide

- Tests in `tools/cftr-model/tests/`: viewer (browser), geometry, research tools, data builder, generated blocks.
- Mutation proofs: 44 deliberate faults (34 in the viewer, 10 in the builder) were put into the code one at a time and the tests run; all 44 were caught. After the off-screen pause was added, eight more viewer faults were written and all 52 were re-run on the final code: 50 caught, 2 survived (two of the three overlapping guards that pause drawing off-screen, whose removal changes nothing observable), none failed to apply. An earlier run's "caught" results were unreliable because a ring-marker test was flaky (it measured 0.875 once in the full suite); the test now counts gap transitions instead of a fraction, and the whole run was repeated. The harness is not in this repository (it lives in the local scratch folder), so this result is a claim from the session, not something a reader can re-run here.
- Data: rebuilt from the fetched files, the shipped data module is byte-identical; 4,609 placed residues were compared with the UniProt sequence and 3 differ, all declared by their files.

## Slips caught on the way (kept so they are not repeated)

| Slip | Caught by | Fixed |
| --- | --- | --- |
| "209 entries at 1,293 positions": 13 entries span several residues and inflated the count | Command output against the entry | Single-residue entries only: 196 at 165 positions, 13 spanning |
| 8EIQ gaps "between 5 and 209 residues" | Recount from the file | "16 to 209" |
| "Every model places" the fit residues | Re-read of the fit rule | "6MSM and that model both place" |
| Provenance wording said the numbers come "by the page's own script" | Review | Numbers come from the build tool; wording corrected, Evidence-page teaser deferred |
| A research sub-agent read `data.rcsb.org/rest`, which robots.txt disallows | Review of the sub-agent's transcript | Logged as a conduct slip; the numbers in the note come from the catalogued downloads, not from that read |
| The page told readers to run `tools/cftr/build_traces.py` and `tests/test_cftr_superpose.py`, paths that no longer exist after the package moved here | A security-review read of the site | The generator now says `tools/cftr-model/build_traces.py` and `tools/cftr-model/tests/test_data.py` |
| Site text said 6MSM shows "the normal protein"; its file records an engineered E1371Q change | Blind review, confirmed in the file | Reworded on the site: "in two of the structures" |

## Open

- A person with CF or a carer has not read the page (T-0081); the draft label stays until they have.
- The standalone page here and the site page share their generated blocks and scripts, but the page text around them is kept by hand in both places. `tools/cftr-model/web/index.html` was updated with the site's reviewed wording on 2026-10-05; changing one means changing the other.
- Layout shift when the viewer replaces its still picture, and the viewer drawing continuously on the home page, were raised by the SEO review from reading the code, not from measurement. Nothing was measured.
- Done in the same release: ledger entries for the slips above (E84 to E92), the README rows and the pointer from `cystic-fibrosis.md`.
