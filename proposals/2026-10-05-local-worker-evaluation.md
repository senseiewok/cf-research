# Evaluation: the local worker on this week's delegations, and what changed in the loop

Status: evaluation, not a policy change. Written 2026-10-05. Everything below is counted from saved attempt logs (`usage.jsonl`, `reply-N.txt`, verifier output); nothing is estimated. The sample is small (ten delegations, one model, one machine), so it describes what happened and does not rank models.

The worker is Qwen3.8 27B (64k context) on Ollama. The loop is `delegate.ps1` in the control repo: the controlling agent writes a self-contained task packet and a verifier first, proves the verifier on a known-bad and a known-good input, sends the packet, and reads every accepted output. A model never marks its own work done.

## What was delegated and what happened

| Task | Kind | Attempts and mode | Outcome | What reading the output found |
| --- | --- | --- | --- | --- |
| `fetch_sources.py --verify` honours `--only` (T-0054) | Python edit, 16 unit tests as the verifier | 1, fast, 17 s | Accepted | Tidied by the controller; behaviour matched the spec |
| Catalog field checker (T-0026) | Python script, fixtures | 1, fast, 12 s | Accepted | Tidied by the controller |
| Quote-completeness check (T-0028) | Python, sentence-boundary arithmetic | 3 rounds of 2 fast and 1 thinking attempt | Accepted on round 3's thinking attempt (15,952 tokens, 105 s) | Early failures were off-by-one offsets. Round 2's thinking attempt hit the 16,384-token cap, with a complete answer inside its thinking text; extracted by hand, it passed the verifier |
| `check-lab-files.ps1` (T-0056) | PowerShell, five rules | 3 rounds, 8 calls | Not accepted; the controller wrote it | My packet put a statement before `param()`; then `$script:` in nested functions; then output mixed with return values |
| Contact validator and terminal-text cleaner (T-0052 A) | Two small Python functions | 1, fast, 12 s, 902 tokens | Accepted | The cleaner's CSI pattern missed `:<=>` parameter bytes; a case was added and the controller rewrote it |
| robots.txt group selection (T-0052 B) | One Python function, 26 cases | 1, fast, 6 s, 560 tokens | Accepted | Logic matched an independent reference on 320,000 random robots files |
| `setup.sh`, three PR descriptions | Shell, Markdown | See the earlier evidence pack | Accepted after repair; PR descriptions accepted with invented framing found by reading | |

Not delegated, with the reason: edits to the delegation driver itself and to the board renderer (the worker cannot test either), catalog entries from web pages (untrusted text and a networked controller are needed), and cross-provider changes whose design was still open.

## Patterns that held

1. **Small, precise, test-shaped tasks succeed on the first fast attempt.** The last three delegations took 6 to 17 seconds each. Each had a spec with worked examples and a verifier with dozens of cases.
2. **A passing verifier does not clear the output.** Reading found a real defect in one of the three most recent accepted outputs and wrong framing in the PR descriptions. The defect was a gap the verifier did not cover; it is now covered.
3. **The spec is the first suspect.** Two of the failures traced to the packet: a statement before `param()` and a test case of mine that contradicted the standard. In the second case the controller's own reference implementation failed the verifier before delegation, which is the point of writing one.
4. **Thinking mode is a budget problem before it is a quality problem.** Thinking attempts spent the whole cap on reasoning and returned no code, while in one case the thinking text already held a complete answer that passed the verifier. The loop now keeps the text and tries its complete answers against the verifier; the one pass that followed came from a thinking attempt that finished under the higher cap.

## Loop changes made from these logs, and their evidence

| Change | Evidence it works | Evidence it does not yet have |
| --- | --- | --- |
| `-ThinkMode` off for model-only runs, explicit `-MaxOutputTokens` up to 16384 | The cap error now names the cap; tested in `test-delegate.ps1` | No measured gain in accepted attempts |
| Salvage complete answers from capped thinking text | A complete answer extracted by hand from the capped round-2 thinking text passed the verifier, which prompted the feature; the code path is unit-tested | Never yet seen to rescue a live run |
| Best-so-far anchor and an identical-resubmission check | Unit-tested; the loop shows the best attempt after a worse rewrite | No second case where it mattered |
| Feedback limited to 12 lines | Prevents the failure list crowding the packet | Not compared with unbounded feedback |

## Too small to conclude

- Whether thinking mode is worth its cost on tasks like these; one salvage is not a rate.
- Whether the local worker is better or worse than a cloud model on the same packets; no paired run was made.
- Whether a three-model council adds anything; an earlier four-task pilot found no gain over a direct fast attempt, and nothing here changes that.

## Proposals that remain proposals

- An AST-level lint for PowerShell candidates (statements before `param()`, `$script:` in nested functions, output mixed with return values) run before the verifier. These three caused the T-0056 failures.
- A two-stage shape for hard tasks: think and write a design, then implement against it, so a capped thinking text is not the only copy of the idea.
- A packet lint that reads the spec as the language would, to catch a wrong spec before a worker obeys it.
- Three experiments from the frontier review of these logs: richer feedback (expected versus actual), worked examples in every packet, and a hand-trace comment requirement. Only the first step of the first (keep the thinking text) is done and confirmed.

## Later the same day: four more delegations, a worker review, and what they changed

| Task | Kind | Attempts and mode | Outcome | What reading found |
| --- | --- | --- | --- | --- |
| Manual-downloads page generator | Python, verifier of fixture and command-line checks | 1, fast, 15 s | Accepted | A dead helper function, removed; a mutation caught it |
| ASCII canvas library | Python, verifier of hand-computed tiny cases | 2 fast attempts, both failing one check | The worker's output was right; my expected value was wrong | A mislabelled algorithm and a duplicated table, tidied |
| ASCII rendering checker | Python, verifier of small Markdown fixtures | 2 fast, then thinking (72 s, 11,507 tokens) | Accepted on the thinking attempt, after I corrected two mistakes in my spec | Clean against its spec |
| Worker review of the T-0052 diff (470 lines, two samples) | Review with the findings schema | 1 call, 2 samples | Three findings, none true | A claimed regular-expression blow-up ran in 0.02 to 0.26 s on 5 MB adversarial input; a claimed escape echo was neutralised by `repr`; one was naming. The controller's own read found one real gap (an unbounded value in an error message) |

Three times in the last two days a failing verifier was the verifier's fault (a line number, a spec gap, a wrong expectation); the loop now says `VERIFIER SUSPECT` when one single check fails with the same message on two different candidates, and `delegate-batch.ps1` runs a list of delegations unattended. Ollama already queues concurrent requests (the call timeout includes the wait), so no lock was added.
