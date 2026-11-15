# Rare-Genotype & Class-I Clinical Trial Screening Benchmark (2026)

*Track 1 Dossier — Clinical Trials & Rare Class-I Genotypes (`T-0940` / `PRED-0003`)*
*Primary Sources:* `cffpr-adr-2024`, `ccfr-adr-2024`, `acfdr-adr-2025`, `ecfspr-adr-2024`, `clinicaltrials-gov`
*Verification Sidecars:* `rare-genotype-trial-screening-2026-claims.json`, `rare-genotype-trial-screening-2026-evidence.txt`, `rare-genotype-trial-screening-2026-evidence-map.json`, `rare-genotype-trial-screening-2026-allow.txt`

---

## 1. Multi-Registry Modulator Ineligibility Denominators

Small-molecule CFTR modulators have transformed outcomes for individuals carrying responsive protein-producing variants, yet national registry denominators confirm a persistent cohort of people with Cystic Fibrosis who remain ineligible by genotype or age:

| Registry Source (`catalog.yaml` ID) | Reporting Year | Total Cohort / Subgroup | Eligible / Prescribed Metric | Ineligible Denominator |
| --- | --- | --- | --- | --- |
| US CFF Patient Registry (`cffpr-adr-2024`) | 2024 | Total genotyped CF population | 92.3 percent genotype-eligible; 26,277 of 29,314 (89.6 percent) age-and-genotype eligible prescribed | 8.0 percent ineligible by genotype for all FDA-approved modulators; 3,037 (10.4 percent) eligible but unprescribed (down from 3,065 or 10.8 percent in 2023) |
| Canadian CF Registry (`ccfr-adr-2024`) | 2024 | Table 7 modulator census | 3,087 individuals on elexacaftor/tezacaftor/ivacaftor (1,004 children and 2,083 adults) | 7.6 percent not considered currently eligible by age and genotype for any CFTR modulator |
| Australian CF Data Registry (`acfdr-adr-2025`) | 2025 | 2,326 adults and 1,473 children/adolescents with known eligibility | 85.6 percent of adults eligible and prescribed | 5.8 percent of adults and 11.5 percent of children/adolescents (170 paediatric patients) not eligible for a modulator |
| European CF Society Patient Registry (`ecfspr-adr-2024`) | 2024 | Participating national registries | Adopted 2024 EMA eligibility criteria across all countries except Israel, the Russian Federation and Switzerland | Non-Class I mutation requirement under EMA Kaftrio/Alyftrek labels |

*Verified Registry Denominator Claims:*
- In 2024, subtracting 26,277 prescribed individuals from 29,314 eligible individuals confirms that 3,037 age-and-genotype eligible individuals had no modulator prescription reported.
- In 2024, subtracting 1,004 children from the 3,087 total individuals on elexacaftor/tezacaftor/ivacaftor in Canada yields 2,083 adults.
- In the Australian registry in 2025, 5.8 percent of 2,326 adults and 11.5 percent of 1,473 children and adolescents (170 paediatric patients) were not eligible for a modulator.

---

## 2. Deterministic Trial Screening & Word-Boundary Disambiguation (`T-0940`)

In the ClinicalTrials.gov API v2 CF snapshot (`clinicaltrials-gov`), 1,053 studies were identified through the condition query and 940 interventional studies were eligible after deterministic exclusions X1 through X6. When screening those 940 eligible interventional studies for rare-genotype, Class-I nonsense/stop, frameshift, splice, and genotype-agnostic nucleic-acid therapies (`tools/trial_atlas/rare_genotype_trial_screener.py`), naive substring matching produced 28 candidate matches, of which 19 were false positives:

- **Pharmaceutical salt collision (`NCT00846781`):** Study `NCT00846781` evaluated inhalation of denufosol tetrasodium and matched naive substring search only because the chemical name `tetrasodium` contains the three-letter substring `aso`.
- **Modulator-only negative control (`NCT00909532`):** Study `NCT00909532` evaluated ivacaftor in subjects aged 12 years and older with the `G551D` gating mutation and serves as a negative control for Class I genetic therapy screening.

Enforcing strict regular-expression word boundaries (`_keyword_matches`) on short acronyms (`aso`, `mrna`, `trna`, `aav`, `lnp`, `ptc`, `tea`, `teae`, `pk`, `cfu`, `lci`) reduced false positives from 19 to 0 (`28 - 19 = 9` true rare-genotype genetic therapy trials across the 940 studies). From the reduction of false positives from 19 to 0 across 940 studies, we infer that regular expression word-boundary guards are essential when screening clinical trial registries for short genetic therapy acronyms.

### 25-Study Labeled ClinicalTrials.gov Benchmark (`PRED-0003`)

Against the 25-study labeled ClinicalTrials.gov genotype-criteria benchmark (`GENOTYPE_CRITERIA_BENCHMARK_25`, consisting of 20 positive studies across 4 classes and 5 negative/modulator-only controls, `25 - 5 = 20`), `rare_genotype_trial_screener.py --benchmark` achieved **1.0 precision**, **1.0 recall**, and **1.0 accuracy** across all 25 studies.

---

## 3. Representative Class-I, Splice-Modulation, and Genotype-Agnostic Clinical Trials

| NCT ID | Intervention & Modality | Target Population | Registry Results Status (2026 Snapshot) |
| --- | --- | --- | --- |
| `NCT05248230` (AEROW) | Aerosolized `4D-710` AAV gene therapy encoding human `CFTR` Delta-R | Adults with CF who are ineligible for or intolerant of CFTR modulators | Phase 1/2 clinical efficacy results are `"not stated"` in the ClinicalTrials.gov registry record |
| `NCT05712538` | Inhaled `ARCT-032` mRNA formulated in lipid nanoparticles | Healthy adult subjects and adults with CF | Phase 1 safety, tolerability, and pharmacokinetics study |
| `NCT06429176` | `SPL84` antisense oligonucleotide splice modulation | Patients with CF carrying the `3849+10kb` C-to-T splicing mutation | Posted clinical efficacy results are `"not stated"` in the ClinicalTrials.gov registry record |
| `NCT04135495` | `ELX-02` translational readthrough therapy | Subjects with CF carrying at least one `G542X` nonsense allele or genetically similar Class I nonsense mutation | Evaluated safety, pharmacokinetics, and sweat chloride / lung function endpoints |

> **Mandatory ClinicalTrials.gov Disclaimer:** ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote United States Federal Government endorsement.
