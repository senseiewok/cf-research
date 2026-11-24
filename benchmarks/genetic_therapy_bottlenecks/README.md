# `benchmarks/genetic_therapy_bottlenecks/` — Bounded Genetic-Therapy Bottleneck Evidence Map (`Item #5`)

Adjudicated 8-study pilot mapping CF genetic and RNA therapy studies across five translational stages (`S1`–`S5`) so cellular expression or registry listings remain visually and structurally separate from pulmonary clinical benefit.

## Stage Separation & Boundaries

1. **`S1_cell_correction_or_expression`** — Cellular gene/transcript editing, transduction, or CFTR expression in the tested system.
2. **`S2_restored_channel_function`** — Measured channel conductance or organoid functional response.
3. **`S3_in_vivo_airway_delivery`** — Administered in vivo airway vector or RNA delivery.
4. **`S4_durability_or_redosing_persistence`** — Multi-dose schedule or persistence beyond transient acute delivery.
5. **`S5_clinical_benefit_or_lung_function`** — Measured human spirometric or pulmonary endpoint in a clinical study.

Every study header records modality, study design, measurement scale, experimental system, controls, follow-up window, and untested downstream steps. Populated cells require typed provenance anchors:

- **`verbatim_passage`:** Repo-relative file path, LF-normalized `source_file_sha256`, line locator, and verbatim `passage_text`.
- **`snapshot_field_pointer`:** Verified ClinicalTrials.gov `snapshot_nct_id`, snapshot field name, and matching substring in [tools/trial_atlas/snapshot.json](../../tools/trial_atlas/snapshot.json).

**Scope limits:** This map contains no efficacy leaderboard or cross-scale ranking, proposes no wet-lab protocols, labels terminated trials separately without equating termination with failed efficacy (`termination_reason_class: unspecified_in_snapshot`), and infers no clinical benefit from unposted registry entries.

## Run Checks

```bash
python benchmarks/genetic_therapy_bottlenecks/validate_bottleneck_map.py --check
python benchmarks/genetic_therapy_bottlenecks/validate_bottleneck_map.py --check --json
python -m unittest discover -s benchmarks/genetic_therapy_bottlenecks -p 'test_*.py' -v
```
