# `benchmarks/infection_outcomes/` — CFF Infection Research Initiative Trial & Diagnostic Outcomes Benchmark (`T-0944` / `PRED-0009`)

Discovery-practice record (`done.md`, `ledger.json`, `failures.jsonl`, and standalone `independent_audit.py`) accompanying `tools/trial_atlas/infection_outcomes_benchmark.py`.

## Run Checks

```bash
python tools/trial_atlas/infection_outcomes_benchmark.py --check
python benchmarks/infection_outcomes/independent_audit.py --check
python ../cf-lab/.claude/skills/discovery-practice/scripts/check_claim_ledger.py --git benchmarks/infection_outcomes/ledger.json
```
