# Cystic Fibrosis Foundation 2025–2026 Research Priorities: Cataloged Source Architecture & Grounded Evidence Synthesis

**Answer Tier:** `T1` / `T2` (peer-reviewed open-access papers, openFDA SPL labels, CFF Patient Registry Highlights 2025, and `ClinicalTrials.gov` snapshot; `www.cff.org` framing page held strictly `[unverified: manual]` per publisher access rules)  
**Companion Scripts & Benchmarks:** [`tools/trial_atlas/cff_priority_radar.py`](../tools/trial_atlas/cff_priority_radar.py), [`tools/trial_atlas/rare_genotype_trial_screener.py`](../tools/trial_atlas/rare_genotype_trial_screener.py), [`tools/trial_atlas/whole_person_visit_prep.py`](../tools/trial_atlas/whole_person_visit_prep.py), [`tools/variant_profile/super_exon_benchmark.py`](../tools/variant_profile/super_exon_benchmark.py)  
**Verification Sidecars:** [`cff-research-priorities-2026-claims.json`](cff-research-priorities-2026-claims.json), [`cff-research-priorities-2026-evidence.txt`](cff-research-priorities-2026-evidence.txt), [`cff-research-priorities-2026-evidence-map.json`](cff-research-priorities-2026-evidence-map.json), [`cff-research-priorities-2026-allow.txt`](cff-research-priorities-2026-allow.txt)  

---

## 1. Non-Clinical Boundary & Epistemic Source Hierarchy

Independent open-source research, not medical advice. Following our 4-Seat Multi-Family Source Council review (`sources/catalog.yaml`), **CF Foundation research priorities are split cleanly by legal permission and verification tier** so an unverified web landing pointer is never elevated to peer-reviewed clinical evidence:

| Catalog ID | What It Covers | Access & Licence | Epistemic Role (`claim_label`) |
| --- | --- | --- | --- |
| `cff-research-priorities` | [CFF Our Research Approach](https://www.cff.org/research-clinical-trials/our-research-approach) framing pointer (`www.cff.org`) | `access: manual`, `robots: unreadable`, `redistribute: false` | `[unverified]` organizational orientation pointer only (`cff-goals.md` human check; never automated, never placed in `*-claims.json`) |
| `cffpr-highlights-2025` | *2025 CFF Patient Registry Highlights* (`34,526` individuals; `2,253` modulator-ineligible; `33,047` without transplant) | `access: manual`, `redistribute: false` | `[verified]` US Registry orientation (page numbers 1–4 verified from local download) |
| `mcgarry-2025-ijns` | US CFF Newborn Screening Consensus Guideline (`DOI 10.3390/ijns11020024`, `PMC12015897`) | `access: manual`, **`CC BY 4.0`** | `[verified]` 7 consensus recommendations (floating IRT, >=95% ancestral sensitivity, full `CFTR2`/sequencing tier) |
| `kraft-2026-bmcmedresmethodol` | ENVISION Clinical-Trial Consent Framework (`DOI 10.1186/s12874-026-02885-2`, `PMC13721633`) | `access: manual`, **`CC BY 4.0`** | `[verified]` 4-part RIVER consent communication model (`Relationship, Information, Values Exploration, Resolution`) |
| `barilla-2026-omtn` | CFF Therapeutics Lab Intron-1 `TI-1` Partial *CFTR* cDNA Exons `2–27` (`DOI 10.1016/j.omtn.2026.103123`, `PMC13717110`) | `access: manual`, **`CC BY 4.0`** | `[verified]` Preclinical basal stem cell integration (`89/92` or `96.74%` canonical Class I nonsense codons under `NM_000492.4`; `17/92` or `18.48%` for Exons `23–27`) |
| `debray-2026-hepatolcommun` | CFF Hepatobiliary Disease Priorities (`DOI 10.1097/HC9.0000000000001061`, `PMC13668289`) | `access: manual`, **`CC BY-NC-ND 4.0`** | `[verified]` Non-commercial attribution-only priority facts; **zero derivatives, reproduced tables, or translations permitted** |
| Subscription abstracts in `cite-check-ledger.jsonl` | Bihler 2026 (`10.1016/j.jcf.2026.08.010`), Georgiopoulos 2026 (`10.1002/ppul.71472`), Putman 2026 (`10.1016/j.jcf.2026.03.021`), Ren 2026 (`10.1001/jamapediatrics.2026.4593`) | Elsevier / Wiley / AMA metadata verified via Crossref & Scite (`retracted: false`, `contrasting: 0`) | `[T1: PubMed abstract only; full text/tables not ingested pending publisher terms]` (`T-0943`, `T-0945`) |

---

## 2. Five Priority Areas Grounded for Lab Tools & Care-Team Guides

1. **Genetic Therapies & Remaining Non-Modulator Cohort (`TI-1` & Theratyping)**:
   - Under `cffpr-highlights-2025` (page 3), among `33,047` individuals without lung transplants in 2025, `10%` received Alyftrek (`186` of `559` newly eligible starting a modulator after the December 20, 2024 approval), `67%` Trikafta, `3%` Kalydeco, `1%` Orkambi, `1%` Symdeko, `11%` were eligible without prescription, `5%` ineligible by genotype, and `2%` ineligible by age.
   - Preclinical *TI-1* partial *CFTR* cDNA integration (`barilla-2026-omtn`, Exons `2–27`, `codons 19..1480` with Exon 1 `c.52..53` and Exon 2 `c.54` completing codon `18`) covers `89/92` canonical Class I nonsense alleles under RefSeq `NM_000492.4` (`PRED-0007`), explicitly excluding Exon `1` variants (`Q2X`, `S4X`, `S13X`), while Bihler et al. (`10.1016/j.jcf.2026.08.010`) evaluated `417` variants *in vitro* for vanzacaftor/tezacaftor/ivacaftor.
   - Active snapshot interventional genetic-therapy records (`16` interventional trials within `tools/trial_atlas/snapshot.json`, `940` PRISMA-included / `1,053` raw interventional) include `NCT06237335` (`RCT2100` mRNA, `ACTIVE_NOT_RECRUITING`), `NCT05248230` (`4D-710` `CFTRΔR` AAV, `RECRUITING`), `NCT05504837` (`KB407` HSV-1, `RECRUITING`), `NCT06526923` (`SP-101`, `RECRUITING`), `NCT07485543` (`BMD003` mRNA, `RECRUITING`), `NCT06429176` (`SPL84` `3849+10kb C->T` ASO, `RECRUITING`), and `NCT06515002` (`BI 3720931` *Lenticlair 1*, `TERMINATED`; extension `NCT06962852`).

2. **Newborn Screening Equity & Genetic-Trial Consent (`mcgarry-2025-ijns`, `kraft-2026-bmcmedresmethodol`)**:
   - `cffpr-highlights-2025` (page 2) documents increasing diversity among newly diagnosed individuals (`13.0%` Hispanic, `3.7%` Black, `3.6%` multiracial, `2.1%` Asian). `mcgarry-2025-ijns` recommends a floating IRT cutoff, >=`95%` variant-panel sensitivity across all ancestral groups, and a `CFTR` sequencing tier, while `kraft-2026-bmcmedresmethodol` supplies the 4-step RIVER consent framework embedded in `rare_genotype_trial_screener.py`.

3. **Hepatobiliary DILI Monitoring & SPL Boxed-Warning Isolation (`debray-2026-hepatolcommun`, `openfda-label`)**:
   - `TRIKAFTA` (`NDA 212273`) and `ALYFTREK` (`NDA 218730`) carry FDA Boxed Warnings for drug-induced liver injury requiring hepatic lab monitoring, whereas `KALYDECO` (`NDA 203188`), `ORKAMBI` (`NDA 206038`), and `SYMDEKO` (`NDA 210491`) mark Boxed Warning as `not stated` in SPL (`AH-4` enforced inside `whole_person_visit_prep.py`). Debray et al. (`10.1097/HC9.0000000000001061`, `CC BY-NC-ND 4.0`) identify non-invasive liver imaging and biomarker knowledge gaps without reproducing third-party text or tables.

4. **Infection Research Initiative (NTM, Phage, & Non-Sputum Diagnostics)**:
   - `cffpr-highlights-2025` records `22%` *P. aeruginosa*, `12%` MRSA, and `12%` nontuberculous mycobacteria (`median age 26`), motivating snapshot interventional studies `NCT06998043` (`BX004` phage), `NCT07275905` (`AchromoPhage`), and `NCT04294043` (`ABATE` IV gallium for NTM; excluded from PRISMA `940` by rule `X4` within raw `1,053`) alongside external observational studies separated strictly inside `cff_priority_radar.py`.

5. **Whole-Person Care: GI, CFRD, Mental Health, Reproduction, & Telehealth**:
   - Orientation metrics in `cffpr-highlights-2025` include `48%` colonoscopy screening (`2021–2025`; baseline age `40`, age `30` post-transplant), `24%` CFRD (`10+`), `29%` adult overweight, `26%` depression and `29%` anxiety (`12+`), `680` pregnancies, `3,337` people aged `50+`, `11,261` telehealth encounters, and median predicted survival of `66` years (`2021–2025`), paired with `NCT06654752` (`STOP-PEDS`).

*Regulatory Registry Reminder (`AH-1`)*: ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
