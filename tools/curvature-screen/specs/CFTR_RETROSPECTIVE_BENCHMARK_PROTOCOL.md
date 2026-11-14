# Draft Protocol: Retrospective CFTR Benchmark of the Curvature (Wronskian) Readout

> **Draft, for the three-bin burst statistic (S1) and the current-noise statistic (S2).** For the whole-sample variability test, `CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md` governs wherever the two differ (false-alarm level, dead-time rule, minimum bursts, held-out split, comparator, and failure criteria).

Status: **DRAFT — not frozen until committed with a hash before any held-out patch is analysed.**
Scope: tests whether a pre-specified curvature statistic carries information about CFTR gating mechanism on existing recordings. It does **not** test channel health or therapeutic rescue.

## 1. Question and allowed claims
**Question.** On existing patch-clamp recordings, does the pre-specified statistic separate ATP-driven wild-type (WT) gating from a hydrolysis-deficient control, and how does it compare in false-positive rate and sensitivity with established likelihood-based dwell-time analysis?

**Allowed conclusions.** "Under assumptions A1–A5 (§6) the statistic is / is not informative about the number and ordering of kinetic steps in gating." **Not allowed:** statements about healthy channel function, drug efficacy, or rescue.

## 2. Data requested (none analysed yet)
Per patch: raw or idealised trace, sampling rate, analogue filter type/cutoff, holding voltage, temperature, [ATP], PKA state, construct (WT; hydrolysis-deficient e.g. Walker-B E1371Q or equivalent), patch ID, number of channels (or maximum simultaneous openings), recording duration, and any rundown notes. A burst-duration list with the same metadata is an acceptable substitute for the single-channel test.
Candidate sources found so far (none downloaded; licence and file contents unverified): Sheppard-lab deposits at data.bris.ac.uk. Csanády-lab raw dwell data not located in an open repository. If neither is available the study needs a data-sharing agreement with a CFTR electrophysiology group.

## 3. Pre-registration rules (all fixed before held-out analysis)
1. **Unit of analysis = patch.** Patch IDs are hashed (SHA-256 of ID + a published salt); patches with hash mod 2 = 0 form the **development set**, the rest the **held-out set**. Only the development set may be inspected.
2. **Event processing.** Idealisation by half-amplitude threshold on the digitised trace; dead time = 2 × filter rise time; critical closed time `t_crit` chosen on the development set by the Colquhoun–Sigworth method and then frozen; burst = openings separated by closures shorter than `t_crit`. Single-channel patches only (max simultaneous openings = 1 over the full record, with a binomial test that a second channel is absent).
3. **Statistic S1 (single channel).** `D2 = ln n1 − 2 ln n2 + ln n3` on disjoint equal bins of burst duration; bin width `w` chosen on the development set from {25, 50, 100} ms and frozen. Also report the nested-survival version at `(w, 2w, 3w)` as a separate, named observable.
4. **Statistic S2 (multichannel).** Autocovariance Wronskian `W_C` on a frozen lag triple, `K = 16` disjoint time blocks per patch, lags ≥ 5 × filter time constant, drift removal by the 4-lag adjacent-difference variant as a pre-declared secondary.
5. **Comparator C1.** Maximum-likelihood dwell-time fit with exponential vs two-step (hypoexponential) models, LRT with parametric-bootstrap critical values; and HMM-based idealisation-free likelihood where available (e.g. QuB/MIL-type analysis).
6. **Decision rule.** A patch is "flagged" if the one-sided statistic exceeds the threshold giving 5 % type-I error on **simulated** exponential nulls matched to that patch's N, filter and dead time. Thresholds are computed on the development set and frozen.

## 4. Outcomes and analysis
- **Primary:** fraction of held-out WT patches flagged (sensitivity) and fraction of held-out control patches flagged (false-positive rate), with patch-level exact binomial CIs; AUC for S1, S2, C1.
- **Comparison:** difference in sensitivity at matched false-positive rate vs C1; McNemar test on paired patches.
- **Dependence:** burst-duration serial correlation (lag-1 Spearman); block bootstrap by patch and by within-patch time segment; report inflation of SE if any.
- **Filter / dead-time sensitivity:** re-run with dead time ×0.5, ×2; report change in flagged fraction.
- **Power plan:** simulated power (`benchmark/single_channel_sample_size.py`) must show ≥ 80 % power at the chosen bin width for the observed N per patch; otherwise the patch is reported as "uninformative", not "negative".

## 5. What counts as failure (pre-declared)
S1 or S2 is declared not useful if, on held-out patches, sensitivity at 5 % false-positive rate is not higher than 50 % of C1's, or if the false-positive rate on controls exceeds 10 %.

## 6. Assumptions (explicit)
- A1 Stationary gating during the analysed segment (rundown removed or segmented).
- A2 Burst durations approximately independent given the patch.
- A3 Idealisation errors (missed events) are small relative to the first bin or modelled by the dead-time sensitivity.
- A4 A single channel per patch for S1; independent identical channels for S2.
- A5 Control construct is kinetically described by a different step structure (to be written down by a CFTR gating expert before analysis; this is a hypothesis, not a result).

## 7. Deliverables
Code (hash-locked), development-set report, frozen thresholds file, held-out report, and a one-page statement of limits. Everything stored with SHA-256 in the package manifest.
