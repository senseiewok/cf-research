# Antibiotic stewardship in CF: a claim check

A short argument for antibiotic stewardship in cystic fibrosis (CF) was checked against the primary papers. It has five claims. Each one below is set against what the source says, with a short exact quote, the link, the date it was read and the limits.

Research note, not medical advice.

Read on 2026-10-07 through the official APIs listed in `sources/catalog.yaml`: Europe PMC REST (`europe-pmc`), NCBI E-utilities (`ncbi-eutils`) and ClinicalTrials.gov API v2 (`clinicaltrials-gov`). Quotes are exact substrings of the text those APIs returned, with markup removed (for example "ppFEV1" is printed with a subscript in the paper). Labels: **observed** means the source says it; **computed** means arithmetic on numbers the source shows; **our inference** means we drew it, not the source.

## Claims table

| # | Claim as given | Verdict | Exact quote(s) | Source and how it was read | Limits |
| --- | --- | --- | --- | --- | --- |
| 1 | STOP2: 277 early responders randomized to 10 or 14 days of IV antibiotics; mean ppFEV1 change 12.8 vs 13.4, difference -0.65 (95% CI -3.3 to 2.0), 10 days not inferior; slower responders no better at 21 than 14 days; 982 randomized | **SUPPORTED** | "Among 982 randomized people, 277 met improvement criteria"; "Mean ppFEV1 change was 12.8 and 13.4"; "a ‒0.65 difference (95% CI [‒3.3 to 2.0])"; "21 days is not superior to 14 days" | Goss et al., *A Randomized Clinical Trial of Antimicrobial Duration for Cystic Fibrosis Pulmonary Exacerbation Treatment*, Am J Respir Crit Care Med 2021, vol 204, pages 1295-1305. https://doi.org/10.1164/rccm.202102-0461oc (PMC8786075, PMID 34469706). Abstract read via Europe PMC REST, 2026-10-07. Enrollment cross-checked on https://clinicaltrials.gov/study/NCT02781610 via API v2, 2026-10-07 | Abstract only; full text, safety data and the registry's results section were not read. Adults only. The trial measured lung function two weeks after treatment, not resistance or later outcomes |
| 2 | A review says people with CF are repeatedly exposed to antibiotics and available strategies often do not eradicate pathogens, helping resistance develop | **SUPPORTED** | "repeatedly exposed to antibiotics"; "frequently inadequate to eradicate the involved pathogens"; "facilitate the development of antimicrobial resistance (AMR)" | Perikleous et al., *Antibiotic Resistance in Patients with Cystic Fibrosis: Past, Present, and Future*, Antibiotics (Basel) 2023, vol 12, article 217. https://doi.org/10.3390/antibiotics12020217 (PMC9951886, PMID 36830128). Abstract and full text read via Europe PMC REST, 2026-10-07 | A narrative review, not a study. It gives no number for how much harm resistance causes people with CF in the parts we searched (see below) |
| 3 | Culture-based susceptibility testing cannot capture the lung microenvironment or biofilm, so a *sensitive* result may not predict response (attributed to PMC10970943) | **PARTLY** | From PMC9951886, not PMC10970943: "cannot account for the lung microenvironment"; "such as biofilm formation"; "does not steadily prejudge a poor clinical result". From PMC10970943: "discrepancy between laboratory studies (in vitro)" and "real-world conditions (in vivo)" | The wording matches Perikleous et al. 2023 (PMC9951886, above). Vitiello et al., *The Impact of Antimicrobial Resistance in Cystic Fibrosis*, J Clin Med 2024, vol 13, article 1711. https://doi.org/10.3390/jcm13061711 (PMC10970943, PMID 38541936). Both full texts read via Europe PMC REST, 2026-10-07 | The claim cites the wrong review. Neither review measures how well susceptibility results predict response; both are reviews without new data. PMC9951886 also says early therapy should be "based on the in vitro susceptibility" |
| 4 | Antibiotic therapy can distort the gut microbiota in CF; Australian infants with CF (preventive antibiotics common) had lower gut diversity than US infants; early gut changes linked to growth | **PARTLY** | "AUS-CF infants had lower stool alpha diversity"; 98% of Australian infants "routinely prescribed antibiotic prophylaxis"; "none of the US-CF infants received prophylaxis"; "Gut microbial diversity was not linked to growth."; "These are associations and causality cannot be determined." | Deschamp et al., *The association between gut microbiome and growth in infants with cystic fibrosis*, J Cyst Fibros 2023, vol 22, pages 1010-1016. https://doi.org/10.1016/j.jcf.2023.08.001 (PMC10840679, PMID 37598041). Abstract via Europe PMC REST; full text via NCBI E-utilities (db=pmc) after Europe PMC returned a server error (500) for it; both 2026-10-07 | Observational; country and prophylaxis cannot be separated. The growth part of the claim is contradicted for diversity. That this is the paper behind the news summary is our inference; the news page was not opened |
| 5 | After CFTR modulators, Pseudomonas prevalence has fallen in registries but chronic infection usually persists; researchers call for prospective studies on stopping inhaled antibiotics | **SUPPORTED** (Burgel et al.); Muhlebach piece **UNVERIFIED** | "decreased in CF registries since the introduction"; "clinical observations suggest"; "chronic P. aeruginosa infections usually persist"; "need for prospective studies"; "consequences of stopping inhaled antibiotic therapy" | Burgel et al., *Considerations for the use of inhaled antibiotics for Pseudomonas aeruginosa in people with cystic fibrosis receiving CFTR modulator therapy*, BMJ Open Respir Res 2024, vol 11, e002049. https://doi.org/10.1136/bmjresp-2023-002049 (PMC11086488, PMID 38702073). Abstract read via Europe PMC REST, 2026-10-07 | A discussion piece by experts, not a study; "usually persist" rests on "clinical observations", not on numbers shown in the abstract. The DOAJ id given with the claim was not checked (DOAJ is not in the catalog). Full text not read |

## Notes per claim

### 1. STOP2

- **Observed.** The results paper is Goss et al. 2021 in the *American Journal of Respiratory and Critical Care Medicine* (PMC8786075). PMC5742044 is the design paper: Heltshe et al., *Study design considerations for the Standardized Treatment of Pulmonary Exacerbations 2 (STOP2)*, Contemporary Clinical Trials 2018, which Europe PMC tags as a clinical trial protocol.
- **Observed.** All numbers in the claim match the results abstract: 982 randomized, 277 early responders, 12.8 and 13.4, ‒0.65 with 95% CI ‒3.3 to 2.0. The paper says this result is "excluding the predefined noninferiority margin", and concludes 10 days "is not inferior to 14 days".
- **Observed.** For slower responders the abstract reports 3.3 and 3.4 mean ppFEV1 changes for 21 and 14 days, a difference of ‒0.10 (‒1.3 to 1.1).
- **Observed.** ClinicalTrials.gov lists NCT02781610 as randomized, interventional, completed, with an actual enrollment of 982.
- **Computed.** 277 + 705 = 982, matching the abstract's split ("the remaining 705").
- Wording: the claim is accurate. Its *no better* can be stated more precisely with the paper's own term, "not superior".

### 2. Repeated exposure and resistance

- **Observed.** The abstract of PMC9951886 says people with CF are "repeatedly exposed to antibiotics" and that current strategies are "frequently inadequate to eradicate the involved pathogens" and "facilitate the development of antimicrobial resistance (AMR)".
- **Observed, numbers on harm.** In a keyword search of the full text (mortality, outcome, predict, susceptibility, culture), the review links resistant organisms with worse disease in words (for example, it says MRSA brings "higher mortality") but gives no figure for harm caused by resistance in people with CF. Not stated.
- Wording: accurate. The source's "frequently inadequate" is close to the claim's *often do not eradicate*.

### 3. Susceptibility testing and response

- **Observed.** The claim's wording comes from the PMC9951886 abstract: susceptibility tests "cannot account for the lung microenvironment" and adaptive mechanisms "such as biofilm formation". Its full text adds that giving a drug the organism resists in the lab "does not steadily prejudge a poor clinical result" (and vice versa), and it notes "discordance between in vitro and in vivo efficacy".
- **Observed.** PMC10970943 says something related but weaker: it names a "discrepancy between laboratory studies (in vitro)" and "real-world conditions (in vivo)", and says "biofilms create a microenvironment" that protects bacteria. We did not find a sentence in it about culture-based susceptibility testing failing to capture the lung.
- **Not stated.** Neither review gives a measure (for example, the share of *sensitive* results followed by a good response) of how well susceptibility predicts response in CF.
- **Suggested wording:**

  > A 2023 review (Perikleous et al., PMC9951886) says culture-based susceptibility tests cannot account for the lung microenvironment or adaptive mechanisms such as biofilms, and that lab resistance does not reliably predict a poor clinical result, or vice versa. It gives no measure of how well testing predicts response.

### 4. Gut microbiota in infants

- **Observed.** Study design: "prospective, observational study" of infants with CF followed for 12 months at four sites in the US and Australia; 78 infants with CF enrolled, a subset (CF N = 40, non-CF disease controls N = 10) gave stool for microbiome analysis.
- **Observed.** Australian infants had lower stool diversity (p < 0.001); 98% of them got daily prophylaxis (amoxicillin-clavulanate) and none of the US infants did. The authors write that prophylaxis "had marked effects on the gut microbiome".
- **Observed, growth.** "Gut microbial diversity was not linked to growth." Australian infants "had better weight gain over time" and higher mean weight-for-age z-scores (p = 0.02). Malnutrition "was associated with depleted Lactococcus"; prophylaxis and malnutrition were linked with predicted lower activity of short-chain fatty acid pathways (predicted from 16S data, not measured).
- **Observed, causation.** "These are associations and causality cannot be determined." The authors also say it was hard to tell whether a Proteobacteria difference came from "antibiotic prophylaxis or continental differences", because all but one Australian infant got prophylaxis.
- **Suggested wording:**

  > In an observational study of infants with CF (Deschamp et al. 2023), Australian infants, almost all on preventive antibiotics, had lower gut microbial diversity than US infants, who received none. Diversity was not linked to growth; the Australian infants gained more weight. The study shows associations, not cause, and cannot separate antibiotics from country.

### 5. Pseudomonas after CFTR modulators

- **Observed.** Burgel et al. 2024 says *P. aeruginosa* prevalence "decreased in CF registries since the introduction" of CFTR modulators, but "clinical observations suggest" that "chronic P. aeruginosa infections usually persist". It highlights a "need for prospective studies" on the "consequences of stopping inhaled antibiotic therapy".
- **Observed.** Europe PMC tags it as a discussion article. It is an expert view, not a trial or registry analysis.
- **Unverified: the Muhlebach piece.** The closest match found is Muhlebach et al., *Changes in factors associated with inhaled antibiotic prescriptions for people with cystic fibrosis over time in the U.S.*, J Cyst Fibros, vol 24, pages 98-104, https://doi.org/10.1016/j.jcf.2024.09.017 (PMID 39389810), abstract read via Europe PMC REST on 2026-10-07. Europe PMC lists its year as 2025. It is a US registry cohort for 2011 to 2019 and finds the share of people with "chronic and intermittent Pa decreased", with changes "even prior to triple-modulators". It does not, in its abstract, say chronic infection persists on modulators or call for stopping studies. Whether it is the piece the claim meant is not confirmed.
- **Suggested wording:**

  > A 2024 expert discussion paper (Burgel et al., BMJ Open Respiratory Research) says Pseudomonas prevalence has fallen in CF registries since CFTR modulators arrived, but that clinical observations suggest chronic infection usually persists, and it calls for prospective studies of stopping inhaled antibiotics.

## What is not covered

Only abstracts were read for STOP2, Burgel et al. and Muhlebach et al.; their full texts, supplements and the STOP2 registry results section were not. The two reviews and the infant gut paper were searched by keyword, not read line by line, so a sentence giving numbers on harm could have been missed. The news summary at cysticfibrosisnewstoday.com and the DOAJ record were not opened, because neither site is in the source catalog; the link between the news item and Deschamp et al. is our inference. No second model checked these verdicts, and no systematic search was done for other studies that agree or disagree. Nothing here assesses whether any treatment approach is right for anyone.

Research note, not medical advice.
