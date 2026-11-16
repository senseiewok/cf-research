# `benchmarks/theratype/` — Open CFTR Theratype Benchmark (`T-0925`)

Versioned, fixed-split (`2026-11-15.v1`), leakage-audited benchmark evaluating whether open variant-effect and genomic foundation models predict CFTR modulator responsiveness across missense, Class I nonsense/frameshift, and splicing variants—without any patient or PHI data.

## Run Checks

```bash
python -m unittest discover -s benchmarks/theratype -v
python benchmarks/theratype/cftr_theratype_benchmark.py --check --falsify
```
