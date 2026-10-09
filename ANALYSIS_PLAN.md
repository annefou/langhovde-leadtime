# Analysis plan: lead time of basal signals before glacier acceleration (exploratory, n = 1)

Written and committed on 2026-10-09, **before any lead-time, onset, correlation or model
result was computed**. Only data inspection had been done: sampling intervals, data gaps
and variable units of the inputs listed in § 0. The prior chain's published numbers
(`nanopubs/imported/CHAIN_SUMMARY.md`, peak timing at centred smoothing) were known
before writing; nothing else about timing was. The git history of this file is the
record of what was decided before results were seen. Later changes go in § Amendments,
dated, with a reason, and never overwrite the text above them.

Context, scope and ground rules: `DISCOVERY.md` (private; the question is an external
researcher's and is not attributed here).

**Standing caveat, to appear in every output.** The proposed test is held-out prediction
of acceleration across events. This record has **one** acceleration event with basal
pressure (Period II, 2–6 Jan 2022). The predictive test cannot be run on it. Everything
below is exploratory and descriptive, labelled **n = 1**.

## 0. Inputs (Step 0)

No re-derivation of anything the prior chain already established.

- **Code:** `annefou/langhovde-meltwater-replication` v1.0.2, Zenodo
  [doi:10.5281/zenodo.23257925](https://doi.org/10.5281/zenodo.23257925), zip MD5
  `d9ee8c8326eb5322b00f80abc718ff64`. Its `notebooks/gnss.py` (`process_track`,
  `local_regression`, `GnssSettings`, `AUTHORS_FIG3`) is imported from the extracted
  archive, not copied or rewritten.
- **Data:** the Mendeley deposit [doi:10.17632/8wvtxg53ry.1](https://doi.org/10.17632/8wvtxg53ry.1),
  each file checked against the SHA-256 published by the Mendeley API, parsed exactly as
  in v1.0.2 `02_data_clean.py` (GNSS `LG05.dat` = GNSS1, `LG04.dat` = GNSS2;
  `pressure_eliminate1.dat` = BH2201; `aws211219_220206.dat`).
- **Deviation from `DISCOVERY.md` Step 0:** the Zenodo record does not contain
  `results/gnss_tracks.nc` or `results/water_levels.nc` (they were not committed).
  They are regenerated here from the deposit with the v1.0.2 code and settings, and
  checked against the v1.0.2 numbers in `results/sensitivity_c3.csv` and the prior
  Outcome (BH2201 peak 51.0 m a.s.l. at 2022-01-03 23:07 UTC; GNSS1 speed peak 2 h
  earlier at h = 12 h). A mismatch stops the analysis.
- **BH2201 water level:** `(p − 9.2) × 1.0197 − 436.1` m a.s.l. (v1.0.2 constants). The
  first 6 h after installation (until 2022-01-01 00:00 UTC) are excluded, as in v1.0.2.
- **All times UTC.** AWS is UTC (v1.0.2 diurnal check); loggers and GNSS as in v1.0.2.
- BH2202 (from 7 Jan) and BH2203 (from 24 Jan) start after Period II and are not used
  for onsets. AWS air pressure is embargoed; no melt model beyond what is listed below.

## 1. Causal (one-sided) processing (Step 1)

**Why.** The v1.0.2 smoother is a centred Gaussian kernel: the value at t uses data up to
about 3 × h after t. A rise can then appear before it happens, and two series with different
shapes are advanced by different amounts. Lead time must be measured on series that use
only the past.

**Grid.** All series are evaluated on one 15 min UTC grid (the GNSS sampling interval),
00/15/30/45 min.

**Bandwidths.** h = 1, 3, 6, 12 h for every series, centred and causal.

**Centred (reference).** v1.0.2 `process_track` with v1.0.2 settings
(`AUTHORS_FIG3` at 12 h, i.e. centred difference at GNSS1 and forward at GNSS2;
`GnssSettings(bandwidth_h=h)` otherwise), but with `step_h = 0.25`. Check: at the full
hours, the 15 min centred track must equal the v1.0.2 hourly track to 1 mm d⁻¹ in speed
and 1 mm in uplift, away from segment ends (± 3h).

**Causal.** The same Gaussian-kernel local linear regression with weights set to zero for
samples later than the evaluation time (half-Gaussian, σ = h). The fit at t uses only
samples with time ≤ t. Implemented as a new function next to an unmodified import of
v1.0.2 `local_regression`; a unit test asserts that changing any sample after t does
not change the output at t, and that the centred branch reproduces v1.0.2.

- **Position, uplift:** the causal fit's intercept at t (x, y, z relative to the first
  sample, as in v1.0.2).
- **Speed:** horizontal magnitude of the **backward** difference of causal positions
  over one grid step, `|r(t) − r(t − 15 min)| / 15 min`, in m d⁻¹. Uses only the past.
- **Acceleration (the endpoint):** `a(t) = [v(t) − v(t − Δ)] / Δ`, Δ = max(1 h, h), in
  m d⁻²; same formula on centred speed for the reference.
- **Segments:** v1.0.2 gap rule (`max_gap_h = 12`). A causal value needs ≥ 1 h of the
  current segment behind it; earlier values are left empty.
- **BH2201 level and AWS series:** the same one-sided Gaussian kernel (weighted mean,
  i.e. local constant; the local linear fit is reported as a check) applied to the 1 min
  level and the 10 min AWS series, sampled on the 15 min grid. The raw 1 min level at
  the grid time is also reported, because pressure is far less noisy than GNSS.
- **AWS forcing series:** air temperature; positive temperature `max(T, 0)` (melt
  proxy, °C); rain intensity (mm h⁻¹). Rain onset is also taken from the raw record
  (first non-zero `rain_intensity` after 2022-01-02 00:00), unsmoothed.

**Reported for Step 1.** Centred and causal series side by side (figure), and the time
shift of each series' Period II peak between centred and causal at each h.

## 2. Onset timing in Period II (Step 2)

**Series.** AWS positive temperature, AWS rain intensity, BH2201 level, GNSS1 uplift,
GNSS1 speed, GNSS2 speed, GNSS1 acceleration. (GNSS2 is afloat; its uplift is tidal and
is not an onset series.)

**Fixed windows** (from the paper text and the deposit start dates, not from our data):
- reference (pre-event) window **R = 2022-01-01 00:00 to 2022-01-02 00:00 UTC** (24 h; the
  only pre-event day with pressure). Sensitivity: R′ = 2022-01-01 12:00 to 2022-01-02 00:00.
- search window **S = 2022-01-02 00:00 to 2022-01-06 00:00 UTC**.

**Onset rules.** Each is applied to each series × {centred, causal} × h ∈ {1, 3, 6, 12} ×
{R, R′}:
- **T(k): threshold.** Reference mean μ and SD σ in R. Onset = first time in S at which
  the series exceeds μ + kσ and stays above it for 1 h (5 consecutive grid points).
  k ∈ {2, 3, 5}. Upward exceedance for every series.
- **CP: change point.** Hinge model on [start of R, time of the series' maximum in S]:
  constant up to τ, then linear. τ on the 15 min grid, least squares. Onset = τ̂.
  Interval: all τ with SSE(τ) ≤ SSE_min × (1 + F₀.₉₅(1, n − 3)/(n − 3)). This interval
  ignores autocorrelation and is too narrow; it is reported as such.
- **Peak:** time of the maximum in S (for comparison with the prior chain).

An onset that would fall inside R (series already rising before 2 Jan) is reported as
"before S" and is not moved.

**Timing resolution.** ± 7.5 min for GNSS-based series (15 min sampling), ± 0.5 min for
BH2201, ± 5 min for AWS; combined in quadrature for a difference.

**Lead.** For each basal series B ∈ {BH2201 level, GNSS1 uplift} and endpoint E ∈
{GNSS1 speed, GNSS1 acceleration, GNSS2 speed}: lead = onset(E) − onset(B), in hours.
**Positive lead = the basal series changes first.** The same is reported for the AWS
series against E and against B.

**Decision rule (on causal series; centred reported alongside).** Over all rule × h × R
combinations for a pair (B, E):
- **"basal leads"** if lead > 1 h in more than half the combinations where both onsets exist;
- **"simultaneous"** if |lead| ≤ 1 h in more than half;
- **"basal lags"** if lead < −1 h in more than half;
- **"indeterminate"** otherwise, or if fewer than half the combinations give both onsets.

The 1 h tolerance is four GNSS samples, set from the sampling, not from results. The
label is always written with "(n = 1 event, one site)". It describes this event; it is
not evidence of a general lead.

**Period I (context only).** The same onset rules for AWS, GNSS1 uplift and GNSS1/GNSS2
speed, with R = 2021-12-20 00:00 to 2021-12-21 00:00 and S = 2021-12-21 00:00 to
2021-12-26 00:00. No pressure exists. GNSS2 starts 23 Dec (no R; T(k) not applicable,
CP only if ≥ 12 h precede its maximum). GNSS1 has a gap from 2021-12-23 22:45 to
2021-12-25 12:20; an onset or maximum adjacent to it is flagged.

## 3. Lead–lag structure (Step 3)

- **Pairs:** BH2201 level and GNSS1 uplift (X) against GNSS1 speed, GNSS2 speed and
  GNSS1 acceleration (Y); also AWS positive temperature and rain intensity against X and Y.
- **Series:** causal (primary) and centred, at each h; 15 min grid; each series linearly
  detrended within the window.
- **Lags:** −24 h to +24 h in 15 min steps. Pearson r(ℓ) = corr(X(t), Y(t + ℓ)).
  **Positive ℓ = X leads Y.** Reported: ℓ* = argmax r and r(ℓ*), and r(0).
- **Windows:**
  - W1, Period II: 2022-01-02 00:00 – 2022-01-06 00:00;
  - W2, common record GNSS1: 2022-01-01 00:00 – 2022-01-12 08:45 (GNSS1 gap starts then);
  - W3, common record GNSS2: 2022-01-01 00:00 – 2022-01-23 13:00;
  - W4, without the event: 2022-01-06 12:00 – 2022-01-12 08:45 (GNSS1), to show how much
    the one event drives W2.
- **Stability:** ℓ* across h and across windows, tabulated. No significance test is
  claimed: the series are autocorrelated and one event dominates.

## 4. A descriptive version of the three nested models (Step 4)

This is **not** held-out prediction across events. It shows what the proposed test
would look like and which terms carry information here.

- **Target:** change in GNSS1 speed over horizon H, Δv_H(t) = v(t + H) − v(t), with
  H ∈ {1, 3, 6 h}. Also GNSS2 speed (M3 then uses GNSS1 uplift and BH2201).
- **Predictors** at t, all causal, at lags 0, 3 and 6 h, at **h = 3 h primary** (1, 6 and
  12 h reported):
  - **M1:** speed history v(t), v(t − 3 h), v(t − 6 h);
  - **M2:** M1 + positive temperature and rain intensity;
  - **M3:** M2 + BH2201 level and GNSS1 uplift.
  Intercept in each; OLS (statsmodels); hourly sampling of t to limit autocorrelation.
  M1 has 4 parameters, M2 10, M3 16.
- **4a, in-sample** on W2 (GNSS1) and W3 (GNSS2): R², adjusted R², AIC, BIC per model.
- **4b, rolling origin:** expanding window; the first forecast origin is the first hour
  at which the training set has ≥ 80 points (5 × M3's parameters); refit every hour;
  forecast Δv_H. RMSE per model, and skill = 1 − RMSE(M) / RMSE(M1). Reported for the
  whole test span and for test origins inside Period II only. It is stated in the output
  whether Period II's onset falls before the first forecast origin (in which case the
  rolling-origin test cannot forecast the event at all).
- **4c, event hold-out (n = 1):** train on the record after the event — 2022-01-06 12:00
  to 2022-01-12 08:45 (GNSS1) and to 2022-01-23 13:00 (GNSS2) — and forecast
  2022-01-01 00:00 to 2022-01-06 12:00. RMSE and skill as in 4b, overall and over the
  onset part of the event (from the CP onset of the endpoint minus 12 h to its peak).
  The training period contains no rain event, so models M2/M3 may be extrapolating.
- **Decision rule.** M3 "adds information here" if its 4c skill relative to M2 exceeds
  0.10 at H = 3 h for GNSS1 at h = 3 h **and** has the same sign at the other two
  horizons. That label is about this event and these data only.

## 5. The knowledge gap (Step 5, main expected output)

A written section, not a computation: how many events, which sensors, sampling and
continuity a proper held-out test needs, sized from what Steps 2–4 show (e.g. the onset
separation that would need to be resolved, and the sensor gaps that hit event onsets
here), and the candidate data sources listed in `DISCOVERY.md`.

## Outputs

- `results/series_15min.nc` — all centred and causal series (NetCDF).
- `results/onsets.csv`, `results/leads.csv`, `results/lead_labels.csv` (Step 2).
- `results/xcorr.csv` (Step 3).
- `results/nested_models.csv` (Step 4).
- `figures/main_result.png` — Period II, causal series and onsets.
- `SUMMARY_PRIVATE.md` — short summary for Anne (decision point in `DISCOVERY.md`).

Every result is labelled pre-registered or post hoc.

## Amendments

**A1 — 2026-10-09, before any result was computed.** § 1 specified a one-sided
local-constant kernel (weighted mean) for BH2201 and the AWS series, with the local
linear fit as a check. Swapped: **all series use the same Gaussian-kernel local linear
fit** (v1.0.2 `local_regression` centred; its one-sided version causal), and the local
constant fit is the check. Reason: a one-sided local-constant kernel lags a linear ramp
by h·√(2/π) ≈ 0.8 h, whereas the local linear fit does not. Using it for pressure but
not for GNSS would have delayed pressure relative to speed by construction, biasing the
comparison against the basal series.
