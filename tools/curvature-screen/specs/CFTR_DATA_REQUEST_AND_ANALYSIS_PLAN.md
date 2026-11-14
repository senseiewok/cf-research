# CFTR Data Request and Analysis Plan

**Status:** fixed on `2026-11-13`, before any CFTR recording has been seen, except that four pieces of code are not yet written: the comparator of step 6, the control power simulation of step 7, the `0.36 / f_c` default of step 2, and the one-channel check of step 1. This plan covers the whole-sample variability test. Where it differs from the earlier draft protocol for the three-bin and current-noise statistics (`CFTR_RETROSPECTIVE_BENCHMARK_PROTOCOL.md`), this plan governs: one-sided `α = 0.01` instead of `5%`; the dead-time rule in step 2; the burst minimums in step 1; no held-out split (step 8); the comparator (step 6); and the failure criteria (section 5). Changes are listed in section 7. File names in this plan refer to the accompanying code bundle.

## 1. Question

Can the burst durations of wild-type CFTR be a positive mixture of exponentials? For one channel at ideal resolution, a burst that is a stay in a set of states of a reversible gating scheme has such a distribution. The check is one-way: a rejection rules the mixture out; a non-rejection shows nothing.

## 2. Data requested

| Item | Requested | Better still |
| :--- | :--- | :--- |
| Constructs | Wild-type CFTR | Plus one hydrolysis-deficient or non-hydrolytic control (for example E1371Q, D1370N, or K1250A) recorded under comparable conditions |
| Unit of data | One row per burst: patch identifier, burst start time, burst duration | Idealized event lists or raw traces, so that bursts can be rebuilt under one stated rule |
| Per patch | Sampling rate, filter type and cut-off, resolution (dead time) imposed, burst criterion (critical closed time), number of active channels | Temperature, [ATP], PKA, voltage, date |
| Amount, wild type | `100` or more bursts per patch from `3` or more patches (patches with `50` to `99` bursts are still used) | `300` or more bursts per patch |
| Amount, control | Whatever exists; long bursts mean few events per patch | As many patches as available |

No patient data are requested or wanted. Commitments to the contributing laboratory: it sees the results first and can correct them; the data are not redistributed and are deleted on request; use is non-commercial; and credit takes whatever form the laboratory prefers, or none.

## 3. Fixed analysis

1. **Inclusion.** Wild-type patches with one active channel and at least `50` bursts after step 2. Control patches with at least `5` bursts. Patches the contributing laboratory flags as unstable are excluded before analysis, and the exclusions are listed.
2. **Dead time.** A single resolution `t_res` is taken as the value the laboratory imposed on openings and shuttings before bursts were built. If none is stated, `t_res = 0.36 / f_c` for the stated `−3 dB` cut-off `f_c` of a Gaussian filter (about twice the usual half-amplitude dead time; a deliberately conservative default), and this is reported. Bursts shorter than `t_res` are discarded and `t_res` is subtracted from the rest. If durations lie on a frame grid, `t_res` is rounded up to the grid and one further frame is discarded, because the first frame above the resolution is only partly observed.
3. **Statistic.** For each patch, the squared coefficient of variation `CV²` (sample variance over squared mean) of the adjusted burst durations.
4. **Null and threshold.** The critical value for each patch is the `1%` quantile of `CV²` for a single exponential of the same sample size and discretization, from `20,000` Monte Carlo draws. One-sided test, `α = 0.01`. The single exponential was the least favourable case among the five two-component mixtures simulated (`N = 100` and `300`); this is not proven in general.
5. **Combining patches.** Each patch's `CV²` is standardized by the mean and standard deviation of its own null. The combined score is `Σ √N_i z_i / √(Σ N_i)`, and its `1%` critical value comes from `20,000` joint Monte Carlo draws using every included patch's actual burst count. Raw durations from different patches are not pooled.
6. **Comparator.** A generalized likelihood-ratio test of the composite null "positive mixture of up to two exponentials" against "two irreversible steps in sequence", on the same bursts, with its critical value set by parametric bootstrap from the fitted null mixture (`2,000` draws per patch, same dead-time rule and discretization). A plain single-exponential-against-two-step likelihood ratio is reported as a secondary, descriptive result, since it can reject for a mixture as well.
7. **Control.** The same steps on the hydrolysis-deficient construct, where no rejection is expected. The test's false-alarm rate is established by simulation, not by the control, which typically has too few bursts to show it; a non-rejection on the control counts as evidence only to the extent of its power, which will be simulated at the observed burst counts and reported (not yet implemented). If no control is supplied, the result is reported as uncontrolled.
8. **No held-out split.** No parameter of this analysis is tuned on the data, so no split is used.

**What is implemented.** `benchmark/whole_sample/run_on_patches.py` implements the dead-time rule of step 2 with `t_res` supplied by the user, steps 3 to 5, the burst minimums of step 1, and the `1×` and `2×` resolution runs of section 5. Not yet implemented: the one-active-channel check (taken from the laboratory's metadata), the `0.36 / f_c` default, the control power simulation, and the comparator of step 6. Those scripts will be added, and their checksums recorded in section 7, before any data are opened.

## 4. What will be reported, whatever the outcome

- Per patch: number of bursts, resolution used, `CV²`, critical value, decision, and the gap between the shortest burst and the stated resolution, with the chance of a gap that large if the resolution is right.
- Combined result per construct at `1×` and `2×` the resolution, on the patches included at both, with the comparator result beside it.
- Every exclusion and every departure from this plan.
- A plain statement of which reading the data support: the check agrees with likelihood analysis; it agrees but needs more data; or it fails.

## 5. Decision rule and what would count as failure

**Decision rule.** A construct counts as rejecting the mixture only if the combined test rejects at both `1×` and `2×` the resolution, on the patches included at both. Rejection at one resolution only is reported as resolution-sensitive and inconclusive: rejecting at `1×` but not `2×` can mean the stated resolution is too small, and the reverse can mean missed or merged events.

**Failure of the check:**

- The combined test rejects for the control construct at `α = 0.01`.
- On wild type the comparator rejects at `α = 0.01` and the combined `CV²` test does not, while simulation at the comparator's fitted two-step rates and the observed burst counts gives the `CV²` test at least `80%` power.

## 6. Limits known in advance

- The statistic is classical (the randomness parameter of single-molecule kinetics). What is being tested is the calibrated procedure, not the idea.
- Power depends on the gating scheme. In simulation at a `1%` false-alarm rate, with two sequential steps of similar length, `100` bursts give about `90%` power; when one step is five times faster than the other, about a third; with one dominant slow step the burst distribution is close to exponential and the check has little power.
- A resolution stated too small biases the check toward rejection. The script reports how far the shortest burst sits above the stated resolution; its flag has a false-flag rate of about `5%` per patch, has little power when bursts are much longer than the resolution, and is not computed on a frame grid. The `2×` run is the main safeguard.
- Bursts are constructed quantities. The simulations so far assume ideal bursts. Misclassification by the critical closed time, superposition of channels, filtering near the resolution limit, serial correlation between bursts, and rundown have not been simulated.
- Burst durations from a published analysis inherit that analysis's burst criterion and missed-event handling.
- Nothing in this plan concerns therapies, patients, or clinical outcomes.

## 7. Changes

- `2026-11-13`, before any data were seen: revised after independent review. Inclusion minimums by construct; dead-time default; `√N` weights; composite-null comparator; failure criteria; held-out split dropped; a resolution-doubling criterion that fired on valid mixtures replaced by the decision rule in section 5.
- `2026-11-13`, later the same day, before any data were seen: revised after a second review. Frame-grid handling added to step 2 (resolution rounded up and one further frame discarded); the `0.5×` resolution run dropped; a sentence describing a burst-construction simulation that is not in the code bundle replaced by a statement that it has not been simulated; the least-favourable-case sentence in step 4 narrowed to what was simulated; the shortest-burst gap added to the report; the list of unwritten code completed; the decision rule separated from the failure criteria.
