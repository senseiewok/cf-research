# tools/curvature-screen

A statistical screening tool for noisy decay curves and for lists of event durations, with its benchmarks on public data and its simulations for CFTR single-channel recordings.

**Not medical advice.** This is research software and a research note. Every CFTR result here is a simulation. No CFTR recording and no patient data have been analyzed.

## What it is

Three readings of a decay curve at equal steps, `y1`, `y2`, `y3`, give `W = y1*y3 - y2^2`. Any mixture of independent exponential decays with positive weights gives `W >= 0`. A significant negative `W` rules that mixture out. It does not say what replaces it, and a non-negative `W` shows nothing.

The folder holds:

- two closed-form estimators of `W` that remove the bias from squaring a noisy reading and cancel an unknown constant background;
- a whole-sample variability check for event lists (burst or dwell durations), which tests the same null hypothesis and, in simulation, has more power than three bins;
- runs on seven public experimental datasets, with control analyses and a false-alarm calibration;
- simulations of how many CFTR bursts, or how much recording time, each test would need under a stated two-step gating model.

## What it does not show

- It is not validated on CFTR recordings. That needs event-level data shared by a laboratory; see the data request below.
- It does not beat likelihood fitting. In simulation under a correctly specified two-step model, a likelihood-ratio test reaches about 100% power at 300 bursts where the three-bin test reaches 29.5% (`BRIEFING.md`, section "Reconciled CFTR patch-clamp simulation scorecard").
- No single ingredient is new. The briefing lists the prior work for each one.
- Several reported results are weak or unsupported: the `RED_DA` onset sign, the XPCS short-lag dome, the ATP-analog values and the dark-tile control. The briefing gives the details.

## Read first

| File | What it is |
| :--- | :--- |
| [`BRIEFING.md`](BRIEFING.md) | The full note: claims, results, limits, prior work. Also as [PDF](Curvature_Sign_Screening_Briefing_V6_2026-11-13.pdf). |
| [`specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md`](specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md) | What data a CFTR test needs, and the analysis fixed before any recording is seen. |
| [`specs/WHOLE_SAMPLE_METHODS_NOTE.md`](specs/WHOLE_SAMPLE_METHODS_NOTE.md) | The whole-sample variability check, step by step. |
| [`specs/CFTR_RETROSPECTIVE_BENCHMARK_PROTOCOL.md`](specs/CFTR_RETROSPECTIVE_BENCHMARK_PROTOCOL.md) | Draft protocol for the three-bin and current-noise statistics. Not frozen. |

## For laboratories with CFTR single-channel recordings

The open question is whether the whole-sample check rejects a passive mixture on real wild-type CFTR burst durations. A table of burst durations per patch is enough, with each patch's resolution (dead time), filter cut-off, burst criterion and number of active channels; raw traces are not needed. The plan in `specs/CFTR_DATA_REQUEST_AND_ANALYSIS_PLAN.md` states the request and the analysis. To get in touch, use the Contact page of [senseiewok.ai](https://senseiewok.ai/contact/).

`benchmark/whole_sample/run_on_patches.py` applies the dead-time rule and steps 3 to 5 of the fixed analysis to such a table, so a laboratory can run it without sending data anywhere.

## Run it

Python 3 with the two packages in `requirements.txt`.

```bash
pip install -r requirements.txt
python benchmark/single_channel_sample_size.py
python benchmark/macroscopic_patch_clamp_benchmark.py
python benchmark/whole_sample/cv2_burst_test.py
```

The scripts are copied byte for byte from the sealed package, because the briefing quotes the checksums of two of them and of the draft protocol. Two scripts therefore still name files that are not part of this release: a docstring in `benchmark/validate_fluorescence_benchmark.py` and comments in `benchmark/controls/gap1_sigma_lambda.py`.

## What is not included

- The Lean sources for the supporting lemmas, and a cross-domain synthetic test suite. Both share files with a separate study that is not released. The briefing names the lemmas it relies on; a reader cannot re-run that formal check from this folder, and no number in the briefing depends on it.
- An earlier findings log and an engine specification, which the briefing supersedes.
- The raw public datasets. `benchmark/out/` holds small derived caches (histograms and correlation tables); the originals are at the links below. The Dataset 1 results (Zenodo 14848839) have no cache here: re-running them needs the raw download, with its folder given in the environment variable `ZENODO_14848839_RAW`.
- The text differs from the sealed version 6 in a few places: the byline, the references to the files listed above, and wording lowered in a last review on 2026-11-14. No number was changed.

## Public data used

Licences were read from each Zenodo record on 2026-11-14. The derived caches in `benchmark/out/` keep the licence of their source.

| Zenodo record | Depositors (first listed) | Licence | Cache here |
| :--- | :--- | :--- | :--- |
| [14848839](https://doi.org/10.5281/zenodo.14848839) | Gäbler, Jain, Then and others | CC BY-SA 4.0 | none |
| [1249497](https://doi.org/10.5281/zenodo.1249497) | Hugel and others | CC BY 4.0 | `hellenkamp2018_zen1249497_histograms.npz` |
| [7656540](https://doi.org/10.5281/zenodo.7656540) | Wetzker | CC BY 4.0 | `zenodo7656540_leica_stellaris8_tissue_flim.npz` |
| [17207776](https://doi.org/10.5281/zenodo.17207776) | Chaney | CC0 1.0 | `zenodo17207776_nsls2_xpcs.npz` |
| [15476468](https://doi.org/10.5281/zenodo.15476468) | Trybek | CC BY 4.0 | `zenodo15476468_bk_channel.npz` |
| [5567245](https://doi.org/10.5281/zenodo.5567245) | Kobayashi, Ueno, Li and others | CC0 1.0 | `zenodo5567245_f1_atpase.npz` |
| [5533877](https://doi.org/10.5281/zenodo.5533877) | Larobina, Pommella, Philippe and others | CC BY 4.0 | `zenodo5533877_mucus_gel.npz` |
| [4574851](https://doi.org/10.5281/zenodo.4574851) | Anderson and others | CC BY 4.0 | `zenodo4574851_grb180805a.npz` |
| [15278967](https://doi.org/10.5281/zenodo.15278967), [15229515](https://doi.org/10.5281/zenodo.15229515) | Li and others; Hinds, Perley, Sollerman and others | CC BY 4.0 | `zenodo15229515_15278967_cosmic_triad.npz` |
| [7215281](https://doi.org/10.5281/zenodo.7215281) | Savchenko and others | CC BY 4.0 | none |

None of these datasets contains CFTR recordings or patient data. Using them implies no endorsement by their depositors.

## Licence

Code in this folder is MIT ([`../LICENSE`](../LICENSE)). The text of the briefing and the three spec documents is CC0-1.0. The briefing and the methods note state it; it applies to the other two spec documents as well. The derived caches keep the licence of the dataset they come from.

`SHA256SUMS.txt` lists a checksum for every file in this folder.
