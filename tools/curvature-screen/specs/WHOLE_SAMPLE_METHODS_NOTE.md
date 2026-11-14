# A Whole-Sample Variability Check for Burst Durations: Methods Note

**Sensei Ewok's CF Lab** (`https://senseiewok.ai`) · `2026-11-13` · Licence: `MIT` (code), `CC0-1.0` (text). Independent, non-commercial work. Not medical or clinical advice.

> **In brief.** This is a statistical check that can rule out one model of burst durations: a positive mixture of exponentials. It uses a classical statistic. It has been run on simulations and on public F₁-ATPase dwell times. It has **not** been run on any CFTR recording. This note describes the procedure, what was tested, and its limits.

## 1. The check

For one channel at ideal resolution, a burst that is a stay in a set of states of a reversible gating scheme has a duration distribution that is a positive mixture of exponentials. Every such mixture has a squared coefficient of variation `CV² ≥ 1`, with equality only for a single exponential. A burst made of two or more irreversible steps in sequence has `CV² < 1`.

The check is one-way. `CV²` significantly below `1` rules the mixture out. `CV² ≥ 1` shows nothing.

The statistic is the classical randomness parameter of single-molecule kinetics ([Svoboda, Mitra & Block 1994](https://doi.org/10.1073/pnas.91.25.11782); [Schnitzer & Block 1995](https://doi.org/10.1101/sqb.1995.060.01.085)). What is offered here is a worked, calibrated procedure, and none of its steps is new on its own.

## 2. Procedure

1. Discard bursts shorter than the resolution `t_res` and subtract `t_res` from the rest. An exponential stays exponential under this rule, and a positive mixture stays a positive mixture. If durations lie on a frame grid, `t_res` is rounded up to the grid and one further frame is discarded.
2. For each patch, compute `CV²` and compare it with the `1%` quantile of `CV²` for a single exponential of the same sample size and discretization (`20,000` Monte Carlo draws). One-sided test, `α = 0.01`.
3. Standardize each patch's `CV²` against its own null and combine patches as `Σ √N_i z_i / √(Σ N_i)`, calibrated by joint Monte Carlo with the actual burst counts. Raw durations from different patches are not pooled.
4. Repeat at twice the resolution on the same patches. The mixture counts as rejected only if the combined test rejects at both resolutions; otherwise the result is reported as resolution-sensitive or as no rejection.

`benchmark/whole_sample/run_on_patches.py` implements these steps for a table of burst durations. The full plan, with inclusion rules, a likelihood-ratio comparator, and failure criteria, is in `CFTR_Data_Request_and_Analysis_Plan_2026-11-13.pdf`.

## 3. What was tested

Power results are for two irreversible steps in sequence, one-sided `α = 0.01`, `20,000` runs unless stated; the dead-time and patch results use `4,000` runs (`benchmark/whole_sample/cv2_burst_test.json`).

**Power depends on how similar the two steps are.**

| Rate ratio (faster / slower) | Population `CV²` | Power, `100` bursts | Power, `300` bursts | Power, `1,000` bursts |
| :--- | :--- | :--- | :--- | :--- |
| `1` | `0.50` | `94%` | `>99.9%` | `>99.9%` |
| `1.5` | `0.52` | `91%` | `99.98%` | `>99.9%` |
| `3` | `0.625` | `61%` | `97%` | `>99.9%` |
| `5` | `0.72` | `32%` | `76%` | `99.68%` |
| `10` | `0.83` | `10%` | `27%` | `73%` |

A single dominant slow step looks almost exponential and is hard to detect.

**False alarms.** On five two-component mixtures (fast-component weight `0.02` to `0.9`, rate ratio `1.2` to `10`), the false-alarm rate was at most `0.9%` (`N = 100` and `300`). The single exponential was the least favourable case among those simulated; this is not proven in general.

**Dead time.** With a dead time of `5 / 20 / 50 ms` and step rates of `12` and `8 s⁻¹`, power at `300` bursts was `99.9%` / `99.2%` / `83.8%`, and the false-alarm rate stayed near `1%`.

**Patches of different speed.** With `10` patches of `100` bursts whose time scales differ (log-normal, `σ = 0.5`), pooling raw durations cut power to `57.3%`; scoring each patch and combining gave `>99.9%` (`4,000` runs).

**Comparison with a three-bin curvature test.** An earlier three-bin version of this check reached `29.5%` power at `300` bursts and `78.3%` at `1,000` under the `12` and `8 s⁻¹` model, where a likelihood-ratio test of the correctly specified model was near `100%` at `300`. The whole-sample check has not been compared with likelihood fitting below `300` bursts.

**Real data (not CFTR).** Public dwell times of bovine mitochondrial F₁-ATPase rotation ([Zenodo 5567245](https://zenodo.org/records/5567245), Kobayashi et al., CC0 1.0), after the dead-time rule:

| Condition | Dwells | `CV²` | Result |
| :--- | :--- | :--- | :--- |
| `3 mM` ATP (saturating) | `1,599` | `0.42` | Mixture rejected (`1%` critical value `0.89`) |
| `1 mM` ATPγS (slowly hydrolysed analog) | `192` | `0.81` | Not rejected (`p = 0.08`) |
| `300 nM` ATP (binding-limited) | `443` | `1.41` | Not rejected; above `1` (post hoc upper-tail `p = 0.0004`), which suggests heterogeneity |

## 4. Limits

- No CFTR recording has been analyzed, and no public event-level CFTR data was found.
- No passive-system control has been run on event lists.
- The simulations assume ideal bursts, one channel per patch, and independent bursts. Misclassification by the critical closed time, superposition of channels, filtering near the resolution limit, serial correlation, and rundown have not been simulated.
- A resolution stated too small biases the check toward rejection. The script reports how far the shortest burst sits above the stated resolution (a flag with about a `5%` false-flag rate per patch, not computed on a frame grid); the run at twice the resolution is the main safeguard.
- A step-finding or burst-building algorithm that misses short events can lower `CV²`.
- The likelihood-ratio comparator in the analysis plan is not yet written.

## 5. Reproducing the results

`python benchmark/whole_sample/cv2_burst_test.py` regenerates the simulation and F₁-ATPase results (`numpy` only). `python benchmark/whole_sample/run_on_patches.py --self-test` runs the analysis script on simulated patches.
