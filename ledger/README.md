# ledger/

Logs of real errors from working sessions and the stated prediction ledger of Sensei Ewok's CF Lab.

## What's inside

| File | What it is | Status |
| :--- | :--- | :--- |
| [`predictions.jsonl`](predictions.jsonl) | One JSON line per stated prediction: id, date, one-sentence prediction, refutation observation, prior art DOIs, standard baseline method, measured public data source (or 'none'), outcome (`supported`, `refuted`, `could_not_test`, `open`), and outside reader check. | Seeded with PRED-0001 (2026-11-14) |
| [`check_predictions.py`](check_predictions.py) | Deterministic validator for `predictions.jsonl` enforcing Discovery Protocol invariants (`docs/lab-focus-2026-11-14.md`). Standard library only. | Working, audited |
| [`test_predictions.py`](test_predictions.py) | Unit tests for schema validation, invariant enforcement, and CLI execution. | 19 tests, 100% pass |
| [`2026-10-04-first-session.md`](2026-10-04-first-session.md) | First error ledger from session logs. | Working, 2026-10-04 |
| [`2026-10-07-second-session.md`](2026-10-07-second-session.md) | Second error ledger from session logs. | Working, 2026-10-07 |
| [`tally.py`](tally.py) | Script to tally error ledger rows by first detector and error class. | Self-check passing |

---

## The Discovery Protocol (`predictions.jsonl`)

The maintainer established the lab's direction on 2026-11-14 (`docs/lab-focus-2026-11-14.md`):

> *A simulation proposes; only a measurement tests. From now on the research sentinel works on one stated prediction at a time.*

### Invariants enforced by `check_predictions.py`:

1. **One line per prediction (`JSONL`)**: Every line is a valid JSON dictionary with no trailing syntax errors.
2. **Required fields**:
   - `id`: Unique identifier formatted as `PRED-NNNN` (e.g. `PRED-0001`).
   - `date`: ISO calendar date string `YYYY-MM-DD`.
   - `prediction`: One sentence stating the prediction clearly.
   - `refutation_observation`: The specific empirical observation that would refute the prediction.
   - `prior_art_dois`: Array of valid DOI strings (`10.xxxx/...`) documenting prior work.
   - `standard_method`: Standard method the proposed idea must match or beat.
   - `measured_public_data`: Admitted source ID from `sources/catalog.yaml`, or `"none"` if no empirical data is publicly available.
   - `outcome`: Strictly one of `supported`, `refuted`, `could_not_test`, or `open`.
   - `checked_by_outside_reader`: Boolean (`true` or `false`).
3. **Single open prediction invariant**: At most ONE prediction may have `outcome == "open"` at any time. A new prediction may only be opened once the previous prediction has been closed by a test, refutation, or filed data request.
4. **Independent review invariant**: An outcome may not be labeled `supported` without `checked_by_outside_reader: true`. Nothing is called a result on a public page until someone outside the lab has checked it.

---

## How to check and test

```bash
# Validate the prediction ledger
python ledger/check_predictions.py

# Run prediction unit tests
python -m unittest discover -s ledger -v

# Run error ledger tally self-check
python ledger/tally.py --self-check
```
