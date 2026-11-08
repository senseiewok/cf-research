# `landscape/` — Research Dossiers and Sentinel Grounding Reports

This directory holds two groups of Markdown notes:
1. **Topical research dossiers and cross-source syntheses** (structural biology, delivery biophysics, trial endpoint evolution, registry denominators, and variant analyses), alongside their companion `*-evidence.txt` and `*-claims.json` verification files.
2. **Dated Chief Researcher Officer (CRO) mission-grounding reports** (`mission-grounding-YYYY-MM-DD*.md`), recording periodic provenance and task-board audits.

## Topical Dossiers

- [`2026-10-24-multi-source-modulator-and-delivery-verification.md`](2026-10-24-multi-source-modulator-and-delivery-verification.md)
- [`alternative-conductance.md`](alternative-conductance.md)
- [`antibiotic-stewardship-cf.md`](antibiotic-stewardship-cf.md)
- [`cf-repositories.md`](cf-repositories.md)
- [`cftr-structures.md`](cftr-structures.md)
- [`delivery-biophysics.md`](delivery-biophysics.md)
- [`discovery-cross-disease-mechanistic-bridge.md`](discovery-cross-disease-mechanistic-bridge.md)
- [`discovery-cross-domain-stem-cell-delivery.md`](discovery-cross-domain-stem-cell-delivery.md)
- [`discovery-cryoem-3d-truncation-model.md`](discovery-cryoem-3d-truncation-model.md)
- [`discovery-dual-aav-intein-packaging.md`](discovery-dual-aav-intein-packaging.md)
- [`discovery-ptc-trimodal-kinetics.md`](discovery-ptc-trimodal-kinetics.md)
- [`discovery-sort-lnp-aerosol-matrix.md`](discovery-sort-lnp-aerosol-matrix.md)
- [`discovery-suppressor-trna-ptc-kinetics.md`](discovery-suppressor-trna-ptc-kinetics.md)
- [`frontier-discovery-2026-10-18.md`](frontier-discovery-2026-10-18.md)
- [`medical-ai-landscape.md`](medical-ai-landscape.md)
- [`modulator-eligibility-denominators-2026.md`](modulator-eligibility-denominators-2026.md)
- [`multi-registry-eligibility-2026.md`](multi-registry-eligibility-2026.md)
- [`phage-coevolution-dynamics.md`](phage-coevolution-dynamics.md)
- [`ptc-readthrough-nmd-biophysics.md`](ptc-readthrough-nmd-biophysics.md)
- [`registry-data-access.md`](registry-data-access.md)
- [`trial-endpoint-evolution-2026.md`](trial-endpoint-evolution-2026.md)
- [`variants/rare-variants-theratypes.md`](variants/rare-variants-theratypes.md)

## Dated Mission-Grounding Sentinel Reports

- [`mission-grounding-2026-10-12.md`](mission-grounding-2026-10-12.md)
- [`mission-grounding-2026-10-13.md`](mission-grounding-2026-10-13.md)
- [`mission-grounding-2026-10-14.md`](mission-grounding-2026-10-14.md)
- [`mission-grounding-2026-10-15.md`](mission-grounding-2026-10-15.md)
- [`mission-grounding-2026-10-16.md`](mission-grounding-2026-10-16.md)
- [`mission-grounding-2026-10-18.md`](mission-grounding-2026-10-18.md)
- [`mission-grounding-2026-10-19.md`](mission-grounding-2026-10-19.md)
- [`mission-grounding-2026-10-20.md`](mission-grounding-2026-10-20.md)
- [`mission-grounding-2026-10-21.md`](mission-grounding-2026-10-21.md)
- [`mission-grounding-2026-10-22.md`](mission-grounding-2026-10-22.md)
- [`mission-grounding-2026-10-23-wave17.md`](mission-grounding-2026-10-23-wave17.md)
- [`mission-grounding-2026-10-23.md`](mission-grounding-2026-10-23.md)
- [`mission-grounding-2026-10-24-wave19.md`](mission-grounding-2026-10-24-wave19.md)
- [`mission-grounding-2026-10-24.md`](mission-grounding-2026-10-24.md)
- [`mission-grounding-2026-10-25-wave21.md`](mission-grounding-2026-10-25-wave21.md)
- [`mission-grounding-2026-10-25-wave22.md`](mission-grounding-2026-10-25-wave22.md)
- [`mission-grounding-2026-10-25-wave24.md`](mission-grounding-2026-10-25-wave24.md)
- [`mission-grounding-2026-10-28-wave25.md`](mission-grounding-2026-10-28-wave25.md)
- [`mission-grounding-2026-10-29.md`](mission-grounding-2026-10-29.md)
- [`mission-grounding-2026-11-02.md`](mission-grounding-2026-11-02.md)
- [`mission-grounding-2026-11-03.md`](mission-grounding-2026-11-03.md)
- [`mission-grounding-2026-11-06.md`](mission-grounding-2026-11-06.md)
- [`mission-grounding-2026-11-07.md`](mission-grounding-2026-11-07.md)
- [`mission-grounding-2026-11-08.md`](mission-grounding-2026-11-08.md)

## Structured Claim & Evidence Verification Sidecars

Several topical syntheses carry machine-checked claim packets (`*-claims.json`), evidence extracts (`*-evidence.txt`), evidence maps (`*-evidence-map.json`), and allowlists (`*-allow.txt`) verified by `tools/claims/check_claims.py` and `tools/claims/check_anti_hallucination.py`:

| Dossier | Verification Sidecar Files |
| --- | --- |
| [`discovery-dual-aav-intein-packaging.md`](discovery-dual-aav-intein-packaging.md) | `discovery-dual-aav-intein-packaging-claims.json`, `discovery-dual-aav-intein-packaging-evidence.txt`, `discovery-dual-aav-intein-packaging-evidence-map.json`, `discovery-dual-aav-intein-packaging-allow.txt` |
| [`modulator-eligibility-denominators-2026.md`](modulator-eligibility-denominators-2026.md) | `modulator-eligibility-denominators-2026.json`, `modulator-eligibility-denominators-2026-claims.json`, `modulator-eligibility-denominators-2026-evidence.txt`, `modulator-eligibility-denominators-2026-evidence-map.json`, `modulator-eligibility-denominators-2026-allow.txt` |
| [`multi-registry-eligibility-2026.md`](multi-registry-eligibility-2026.md) | `multi-registry-eligibility-claims.json`, `multi-registry-eligibility-allow.txt` |
| [`trial-endpoint-evolution-2026.md`](trial-endpoint-evolution-2026.md) | `trial-endpoint-evolution-2026-claims.json`, `trial-endpoint-evolution-2026-evidence.txt`, `trial-endpoint-evolution-2026-evidence-map.json`, `trial-endpoint-evolution-2026-allow.txt` |
