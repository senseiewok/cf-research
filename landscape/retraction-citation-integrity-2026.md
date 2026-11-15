# Scientific Evidence Integrity, Regulatory Label Provenance & Retraction Radar (2026)

*Track 3 Dossier — Scientific Evidence Integrity & Retraction Surveillance (`T-0942`)*
*Primary Sources:* `openfda-label`, `openfda-drugsfda`, `ema-medicines-json`, `cffpr-adr-2024`, `clinicaltrials-gov`
*Verification Sidecars:* `retraction-citation-integrity-2026-claims.json`, `retraction-citation-integrity-2026-evidence.txt`, `retraction-citation-integrity-2026-evidence-map.json`, `retraction-citation-integrity-2026-allow.txt`

---

## 1. Cross-Jurisdictional Regulatory Label Provenance (`openFDA` vs. `EMA`)

Maintaining zero-hallucination clinical boundaries requires tracking exact prescribing-information updates across both the United States Food and Drug Administration (`openfda-label`, `openfda-drugsfda`) and the European Medicines Agency (`ema-medicines-json`):

| Regulatory Authority & Source ID | Product | Effective Milestone | Verified Indication Boundary | Homozygous Class-I Nonsense Status |
| --- | --- | --- | --- | --- |
| US FDA (`openfda-label`) | `TRIKAFTA` (elexacaftor / tezacaftor / ivacaftor) | 03/2026 label update | Patients aged 2 years and older with a clinical diagnosis of CF and at least one variant that is responsive based on clinical and/or in vitro data (Table 6) or results in production of `CFTR` protein | Clinical efficacy data for homozygous Class I nonsense genotypes without `CFTR` protein production are `"not stated"` in the 03/2026 label |
| US FDA (`openfda-label`) | `ALYFTREK` (vanzacaftor / tezacaftor / deutivacaftor) | 03/2026 label update | Patients 6 years of age and older with a clinical diagnosis of CF and at least one variant that is responsive based on clinical and/or in vitro data (Table 5) or results in production of `CFTR` protein | Clinical efficacy data for homozygous Class I nonsense genotypes without `CFTR` protein production are `"not stated"` in the 03/2026 label |
| US FDA (`openfda-drugsfda`) & `cffpr-adr-2024` | `VTD` (`ALYFTREK` initial approval) | December 20, 2024 | Approved for patients aged 6 years and older with at least one `F508del` mutation or another responsive `CFTR` mutation, including 31 additional rare mutations (135 individuals in the Registry became eligible for the first time) | Excludes genotypes lacking responsive or protein-producing `CFTR` variants |
| EMA (`ema-medicines-json`) | `Kaftrio` | 2024–2026 EMA summary | Patients aged 6 years and older and paediatric patients aged 2 to less than 6 years who have at least one non-Class I mutation in the `CFTR` gene | Explicitly requires at least one non-Class I mutation |
| EMA (`ema-medicines-json`) | `Alyftrek` | 2025–2026 EMA summary | People aged 6 years and older who have at least one non-Class I mutation in the `CFTR` gene | Explicitly requires at least one non-Class I mutation |

---

## 2. Deterministic Retraction & Citation Radar Audit (`T-0942`)

To guarantee that no research dossier or structured claim packet relies on a retracted study, withdrawn preprint, or uncataloged citation, `tools/claims/retraction_citation_radar.py` continuously audits `sources/catalog.yaml`, every Markdown dossier in `landscape/` and `proposals/`, and every structured `*-claims.json` and `*-evidence-map.json` sidecar:

1. **Primary Catalog & PRISMA Denominator Audit:** Across `sources/catalog.yaml`, **347 cataloged primary sources** were audited with **0 retracted or withdrawn entries** detected. In the companion ClinicalTrials.gov API v2 CF snapshot (`clinicaltrials-gov`), 1,053 studies were identified and 940 interventional studies were eligible after excluding 113 studies (`1053 - 113 = 940`) under rules X4 (72 studies) and X5 (41 studies, `72 + 41 = 113`).
2. **Adversarial Falsification Suite (`--falsify`):** The radar evaluates 3 deterministic poisoned negative controls:
   - `F1` (`poisoned_catalog_marked_retracted`): Injects a catalog entry marked `retracted: true`.
   - `F2` (`poisoned_external_update_to_retraction_with_injection`): Injects a Crossref `update-to` retraction notice carrying an adversarial prompt-injection string (`<untrusted_data>` escape attempt) and verifies that the notice is both flagged as a retraction and quarantined inside a nonce-tagged `<untrusted_data_NONCE>` boundary.
   - `F3` (`poisoned_uncataloged_claim_citation`): Injects an uncataloged phantom citation ID and verifies that the radar fails closed.
3. **Grounded Synthesis:** All 3 poisoned negative controls failed deterministically as required (3 of 3 falsification controls passed), and external notice titles were quarantined inside nonce-tagged `untrusted_data` boundaries. From the deterministic detection of all 3 poisoned controls across 347 cataloged sources, we infer that coupling catalog provenance checks with nonce-tagged quarantine prevents retracted or adversarial external notices from corrupting clinical synthesis.
