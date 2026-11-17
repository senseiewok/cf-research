# Proposal: Alyftrek (`NDA 218730`) 417-Variant Theratyping Alignment & Intron-1 Partial *CFTR* cDNA (`Exons 2–27`) Exon-Span Coverage Benchmark (`PRED-0007`)

**Status:** proposed, written 2026-11-16  
**Authors:** Sensei Ewok's CF Lab (Carbon–Argon Double-Helix Synthesis)  
**Answer Tier:** `T1` / `T2` (literature metadata verified via PubMed / Europe PMC / Crossref / Scite batch; exon-span benchmark evaluated deterministically on RefSeq `NM_000492.4` and 92 canonical Class I nonsense alleles)

---

## 1. Non-Clinical Dignity & Publisher-Terms Boundary

Independent open-source research, not medical advice. Preclinical assay panels, *in vitro* Fischer Rat Thyroid (FRT) channel recordings, and RefSeq nucleotide-exon arithmetic describe experimental laboratory systems rather than clinical efficacy or personal treatment eligibility.

---

## 2. Two 2026 Primary Literature Anchors & Licensing Separation

1. **Vanzacaftor/Tezacaftor/Ivacaftor (*in vitro* FRT Comparator to FDA-Approved `ALYFTREK` `NDA 218730` [Vanzacaftor/Tezacaftor/Deutivacaftor]) 417-Variant Theratyping — Bihler et al. (Sep 2026)**:
   - **Citation**: Bihler H, Sivachenko A, Bhatt P, et al. *In Vitro Responsiveness of 417 CFTR Variants to Vanzacaftor/Tezacaftor/Ivacaftor*. *Journal of Cystic Fibrosis* (Sep 2026). `DOI 10.1016/j.jcf.2026.08.010`, `PMID 42728175`. Note that the FDA-approved `ALYFTREK` tablet (`NDA 218730`) substitutes deuterated potentiator **deutivacaftor** (`VX-561`) for standard ivacaftor evaluated in this FRT assay panel.
   - **Retraction & Citation Polarity (`scite-batch`)**: `retracted: false`, `editorial_notices: []`, `contrasting: 0`.
   - **Publisher-Terms Boundary**: Published by Elsevier BV (`Journal of Cystic Fibrosis`); only bibliographic metadata and the public abstract are used here. Bulk supplementary table ingestion is **held pending publisher terms verification** (`T-0943` follow-up) and must not be scraped or redistributed without documented permission, whereas public FDA `NDA 218730` Structured Product Label variant lists remain directly auditable.

2. **Targeted Intron-1 (`TI-1`) Insertion of Partial *CFTR* cDNA (`Exons 2–27`) — Barillà et al. (Dec 2026)**:
   - **Citation**: Barillà C, Suzuki S, Park SH, et al. *Targeted insertion of partial CFTR cDNA replicates native CFTR expression with restored channel function in CF airway epithelium*. *Molecular Therapy: Nucleic Acids* (Dec 2026). `DOI 10.1016/j.omtn.2026.103123`, `PMID 42959355`, `PMCID PMC13717110` (`sources/catalog.yaml` id `barilla-2026-omtn`, `CC BY 4.0` gold open access; alongside CFF research priorities pointer `cff-research-priorities`).
   - **Retraction & Citation Polarity (`scite-batch`)**: `retracted: false`, `editorial_notices: []`, `contrasting: 0`.
   - **Mechanism**: Integrates a partial *CFTR* cDNA cassette spanning **exons 2 through 27** (`TI-1`) into *CFTR* intron 1 in primary airway basal stem cells so transcription initiates from native regulatory elements and splices endogenous **exon 1** onto the therapeutic **exons 2–27** sequence.

---

## 3. Deterministic Preclinical Benchmark Implemented (`PRED-0007`)

We implemented [`tools/variant_profile/super_exon_benchmark.py`](../tools/variant_profile/super_exon_benchmark.py) and its unit test suite [`tools/variant_profile/test_super_exon_benchmark.py`](../tools/variant_profile/test_super_exon_benchmark.py) against RefSeq `NM_000492.4` (`CFTR_EXON_SPANS`, 27 exons, 1,480 codons) and the 92 indexed Class I nonsense alleles (`CANONICAL_92_NONSENSE_ALLELES` from Round 122):

| Strategy Regime | Coding Span (`NM_000492.4`) | Covered Nonsense Alleles (of 92) | Explicit Negative Controls Excluded |
| --- | --- | --- | --- |
| **Intron-1 `TI-1` Partial cDNA (`Exons 2–27`)** (`Barillà et al. 2026`) | Exons `2..27` (`codons 19..1480`; native Exon 1 contributes `c.52..53` and donor Exon 2 completes split codon `18` at `c.54`) | **`89 / 92` (`96.74%`)** | `Q2X`, `S4X`, `S13X` in Exon `1` (`codons 1..17` + 2 nt of `18`) |
| **Distal C-Terminal Cassette (`Exons 23–27`)** (comparison window) | Exons `23..27` (`codons 1240..1480`) | **`17 / 92` (`18.48%`)** | `F508del` (Exon `11`), `G542X` (Exon `12`), `R1162X` (Exon `22`), `W1204X` (Exon `22`) |
| **Full Coding Replacement (`Exons 1–27`)** (`NCT06237335` mRNA, `NCT05504837` HSV-1; cf. `NCT05248230` `CFTRΔR`) | Exons `1..27` (`codons 1..1480`) | **`92 / 92` (`100.0%`)** | Promoter/5'-UTR architecture outside coding span not modeled |

*Reminder (`AH-1`)*: `ClinicalTrials.gov` interventional trial IDs cited above (`NCT06237335` RCT2100, `NCT05504837` KB407, `NCT05248230` 4D-710 encoding `CFTRΔR`) are verified in `tools/trial_atlas/snapshot.json`; ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.

---

## 4. Verification & Proposed Next Step (`T-0943`)

- Deterministic verification: `python -m unittest discover -s tools/variant_profile -v`, `python tools/variant_profile/super_exon_benchmark.py --check`, `python tools/variant_profile/super_exon_benchmark.py --falsify`, `python ledger/check_predictions.py`.
- Proposed follow-up (`T-0943`): Request or verify publisher permissions for Bihler et al. (`DOI 10.1016/j.jcf.2026.08.010`) versus public FDA `NDA 218730` SPL variant lists, and benchmark concordance across modulator-responsive vs. *TI-1* exon-1–27 classes once primary variant-level data are admitted to `sources/catalog.yaml`.
