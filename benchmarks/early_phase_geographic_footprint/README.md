# Early-Phase Multi-Country Clinical Trial Site-Network Asymmetry Benchmark (`PRED-0018`)

Independent open-source research benchmark only; not medical advice. Sensei Ewok's CF Lab builds open tools to examine scientific questions; we do not provide clinical advice, interpret individual genotypes, assess personal clinical trial eligibility or site availability, or make treatment promises. Please discuss any care questions with your Cystic Fibrosis (CF) care team. ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.

In plain terms, we asked whether early-stage inhaled genetic and ribonucleic acid (RNA) trials for CF (cystic fibrosis) are much less likely to open hospital sites in two or more countries than early-stage pill or systemic CFTR (cystic fibrosis transmembrane conductance regulator) trials. Measured ClinicalTrials.gov registry records refuted that prediction (`outcome = refuted`). Nearly half of systemic early-stage trials (`26/54`, `48.1%`) and one-third of pure inhaled trials (`5/15`, `33.3%`, or `7/17`, `41.2%` including dual-route regimens) listed clinical sites across at least two sovereign countries (`p = 0.2350`, above the pre-registered `0.01` threshold), showing that rare-genotype inhaled trials already recruit across international site networks despite operating fewer total clinics per study (median `5.0` vs `11.0`).

---

### Technical Deep Dive: Biophysical & Registry Evidence [Tier T2, Internally Checked]

- **Prediction ID & Pre-Registration:** `PRED-0018`, pre-registered in [`done.md`](done.md) (`done_sha256 = b8db87733513b2a3f44dee07b6011a023612b9f084b01fddf4ebc89b83e75ca4`) prior to fetching `protocolSection.contactsLocationsModule` from ClinicalTrials.gov API v2.
- **Locked Cohort (`71` PRISMA-940 Early-Phase Non-Healthy-Volunteer-Only CFTR Trials):**
  - `N_sys = 54` systemic CFTR modulator and systemic premature termination codon (PTC) readthrough trials (`[observed]` from `benchmarks/delivery_route_compartmentalization/results_2026-11-24.json`).
  - `N_inh_pure = 15` pure airway-restricted CFTR genetic, messenger RNA (mRNA), viral-vector, liposomal, and splice-switching antisense oligonucleotide (ASO) trials, and `N_inh_all = 17` total inhaled CFTR trials adding dual-route `NCT05668741` (`VX-522 + ivacaftor`) and `NCT06237335` (`RCT2100 + ivacaftor`).
- **Measured Outcome (`outcome: "refuted"`, zero post-hoc tuning):**
  - **Pure-route primary comparison (`54` vs `15`):** `26/54` (`0.481481`) systemic vs `5/15` (`0.333333`: `NCT02532764`, `NCT02564354`, `NCT06429176`, `NCT06515002`, `NCT06962852`) pure-inhaled (`risk_difference = +0.148148` `[computed]`, one-sided Fisher hypergeometric `p = 10424417/44368472 = 0.23495100 >= 0.01`, two-sided `p = 17117417/44368472 = 0.38580136`, Chan-Zhang unconditional 95% CI `[-0.149733, 0.396110]`, Newcombe-Wilson Method-10 95% CI `[-0.132089, 0.371460]`).
  - **All-17 inhaled comparison (`54` vs `17`):** `26/54` (`0.481481`) vs `7/17` (`0.411765`), `risk_difference = +0.069717`, one-sided Fisher `p = 98997089/239686202 = 0.41302790`, Chan-Zhang 95% CI `[-0.210126, 0.327149]`, Newcombe-Wilson 95% CI `[-0.191690, 0.304623]`.
  - **Complete-case sensitivity (`>= 1` non-withdrawn country, `53` vs `15`):** `26/53` (`0.490566`) vs `5/15` (`0.333333`), one-sided Fisher `p = 3976723/18326108 = 0.21699768`.
  - **Non-gated descriptive jurisdiction & facility footprint (Non-Pooling Invariant):** Median non-withdrawn facilities per trial is `11.0` (`12.0` complete-case, sum `933`) in systemic trials vs `5.0` (sum `134` pure-inhaled, `183` all-inhaled) in inhaled trials. Pure-inhaled site presence is `8/15` (`53.3%`) in the `United States`, `5/15` (`33.3%`) in the `United Kingdom`, `5/15` (`33.3%`) across `EU-27`, `1/15` (`6.7%`) in `Canada`, and `0/15` (`0.0%`; `1/17` all-inhaled) in `Australia`, juxtaposed without cross-registry pooling alongside national registry contexts (`CFFPR 2024`: `8.0%` ineligible by genotype; `CCFR 2024`: `7.6%` by age and genotype; `UKCFR 2024`: `12.4%` never used modulator among `11,381` registered; `ACFDR 2025`: `5.8%` adult and `11.5%` pediatric non-transplant PBS-ineligible).

### Reproducibility Commands

```powershell
python benchmarks/early_phase_geographic_footprint/acquire_locations_fixture.py --help
python benchmarks/early_phase_geographic_footprint/early_phase_geographic_footprint.py --check --falsify
python benchmarks/early_phase_geographic_footprint/independent_audit.py --check
python -m unittest discover -s benchmarks/early_phase_geographic_footprint -v
```
