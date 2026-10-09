# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 03 — Lead-time analysis (Steps 0–4), exploratory, n = 1
#
# Implements `ANALYSIS_PLAN.md` (pre-registered 2026-10-09, amendment A1). **The
# proposed held-out predictive test cannot be run on this record: it has one
# acceleration event with basal pressure (Period II, 2–6 Jan 2022).** Everything here
# describes that one event at one site.
#
# Sign conventions: lead = onset(endpoint) − onset(basal series), positive = basal
# changes first. Cross-correlation lag ℓ positive = X leads Y.

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr

PRIOR = Path("../data/raw/prior_v1.0.2")
sys.path.insert(0, str(PRIOR / "notebooks"))
from gnss import (AUTHORS_FIG3, GnssSettings, _difference, find_segments,  # noqa: E402  (v1.0.2, unmodified)
                  local_regression, process_track)
from leadtime import (backward_rate, onset_changepoint, onset_threshold,  # noqa: E402
                      smooth_on_grid, xcorr)

CLEAN = Path("../data/clean")
RESULTS = Path("../results")
RESULTS.mkdir(exist_ok=True)
BANDWIDTHS = [1, 3, 6, 12]
STEP_H = 0.25
GRID = pd.date_range("2021-12-19 15:00", "2022-02-06 06:00", freq="15min")
PROCS = ["centred", "causal"]
GRID_D = (GRID.values - GRID.values[0]) / np.timedelta64(1, "D")

# %% [markdown]
# ## Step 0 — inputs and regression check against v1.0.2
#
# BH2201 water level with the v1.0.2 constants, first 6 h after installation removed.
# The v1.0.2 hourly GNSS tracks are re-run with v1.0.2 code and settings and compared
# with the archived `results/sensitivity_c3.csv` (Period II peak time and speed, every
# bandwidth). Any mismatch stops the notebook.

# %%
p = xr.open_dataset(CLEAN / "pressure_BH2201.nc")
level = ((p.pressure - 9.2) * 1.0197 + p.z_bed_m).to_series().dropna()["2022-01-01":]
pk = level["2022-01-02":"2022-01-06 23:59"]
print(f"BH2201 Period II peak: {pk.max():.1f} m a.s.l. at {pk.idxmax()}")
assert round(pk.max(), 1) == 51.0 and pk.idxmax() == pd.Timestamp("2022-01-03 23:07")

gnss = {st: xr.open_dataset(CLEAN / f"gnss_{st}.nc") for st in ["GNSS1", "GNSS2"]}


def prior_settings(st: str, h: float, step_h: float = 1.0) -> GnssSettings:
    s = AUTHORS_FIG3[st] if h == 12 else GnssSettings(bandwidth_h=h)
    return GnssSettings(bandwidth_h=s.bandwidth_h, max_gap_h=s.max_gap_h, step_h=step_h, scheme=s.scheme)


def run_track(st: str, h: float, step_h: float) -> xr.Dataset:
    g = gnss[st]
    return process_track(g.time.values, g.lat.values, g.lon.values, g.z.values, prior_settings(st, h, step_h), st)


ref = pd.read_csv(PRIOR / "results" / "sensitivity_c3.csv")
ref = ref[(ref.gap == "observed") & (ref.windows == "paper") & (ref.period == "II")]
hourly = {(st, h): run_track(st, h, 1.0) for st in gnss for h in BANDWIDTHS}
for row in ref.itertuples():
    v = hourly[row.station, row.h_hours].speed.to_series()["2022-01-02":"2022-01-06 23:59"]
    assert abs(v.max() - row.peak_speed) < 1e-9, (row.station, row.h_hours)
    assert abs(v.idxmax() - pd.Timestamp(row.peak_time)) < pd.Timedelta("1s"), (row.station, row.h_hours)
print(f"v1.0.2 Period II peaks reproduced for {len(ref)} station × bandwidth rows")

# %% [markdown]
# ## Step 1 — centred and causal series on one 15 min grid
#
# Centred GNSS: v1.0.2 `process_track` at `step_h = 0.25` is checked against the
# v1.0.2 hourly track below. Its grid starts at each segment's first fix, which is
# off the UTC quarter-hour (GNSS1 from 2021-12-25 12:20:50), so the series used here
# call the same v1.0.2 `local_regression` (all samples, as `process_track` does) and
# `_difference` (v1.0.2 scheme) at the UTC grid times inside each segment. Causal: one-sided local linear fit;
# speed from backward differences. Acceleration a(t) = [v(t) − v(t − Δ)]/Δ, Δ = max(1 h, h).

# %%
def centred_on_grid(t: np.ndarray, v: np.ndarray, st: str, h: float) -> np.ndarray:
    settings = prior_settings(st, h, STEP_H)
    t0 = GRID.values[0]
    td = (t - t0) / np.timedelta64(1, "D")
    gd = (GRID.values - t0) / np.timedelta64(1, "D")
    out = np.full(GRID.size, np.nan)
    for sl in find_segments(t, settings.max_gap_h):
        gi = np.flatnonzero((GRID.values >= t[sl][0]) & (GRID.values <= t[sl][-1]))
        if gi.size:
            out[gi] = local_regression(td, v, gd[gi], h / 24)
    return out


def gnss_series(st: str, h: float) -> dict[str, np.ndarray]:
    c = run_track(st, h, STEP_H)
    g = gnss[st]
    t = g.time.values
    x, y, z = (g[k].values - g[k].values[0] for k in ("x", "y", "z"))
    cx, cy, cz = (smooth_on_grid(t, v, GRID, h, causal=True) for v in (x, y, z))
    px, py = (centred_on_grid(t, v, st, h) for v in (x, y))
    scheme = prior_settings(st, h).scheme
    return {
        f"speed_{st}_centred": np.hypot(_difference(GRID_D, px, scheme), _difference(GRID_D, py, scheme)),
        f"uplift_{st}_centred": centred_on_grid(t, z, st, h),
        f"speed_{st}_causal": np.hypot(backward_rate(cx, STEP_H), backward_rate(cy, STEP_H)),
        f"uplift_{st}_causal": cz,
        "_check": c,
    }


aws = xr.open_dataset(CLEAN / "aws.nc")
forcing = {
    "tpos": (aws.time.values, np.maximum(aws.air_temperature.values, 0.0)),
    "rain": (aws.time.values, aws.rain_intensity.values),
    "level": (level.index.values, level.to_numpy()),
}

series: dict[str, pd.Series] = {}
checks = []
for h in BANDWIDTHS:
    lag = int(max(1, h) / STEP_H)
    for st in gnss:
        out = gnss_series(st, h)
        c = out.pop("_check")
        hr = hourly[st, h]
        common = np.intersect1d(c.time.values, hr.time.values)
        for var, tol in (("speed", 1e-3), ("uplift", 1e-3)):
            a = c[var].sel(time=common).to_series()
            b = hr[var].sel(time=common).to_series()
            seg = hr.segment.sel(time=common).to_series()
            # away from segment ends (± 3h), as pre-registered
            span = seg.groupby(seg).transform(lambda s: s.index.min()), seg.groupby(seg).transform(lambda s: s.index.max())
            inner = (a.index >= span[0] + pd.Timedelta(hours=3 * h)) & (a.index <= span[1] - pd.Timedelta(hours=3 * h))
            diff = float((a - b)[inner].abs().max())
            checks.append({"station": st, "h": h, "var": var, "max_abs_diff": diff, "tolerance": tol, "ok": diff <= tol})
        for k, v in out.items():
            series[f"{k}_h{h}"] = pd.Series(v, index=GRID)
        for proc in PROCS:
            v = series[f"speed_{st}_{proc}_h{h}"].to_numpy()
            series[f"accel_{st}_{proc}_h{h}"] = pd.Series(backward_rate(v, STEP_H, lag), index=GRID)
    for name, (t, v) in forcing.items():
        for proc in PROCS:
            series[f"{name}_{proc}_h{h}"] = pd.Series(
                smooth_on_grid(t, v, GRID, h, causal=proc == "causal"), index=GRID)
        series[f"{name}_causal-const_h{h}"] = pd.Series(
            smooth_on_grid(t, v, GRID, h, causal=True, degree=0), index=GRID)
series["level_raw"] = pd.Series(level.reindex(GRID, method="ffill", tolerance=pd.Timedelta("2min")).to_numpy(), index=GRID)

checks = pd.DataFrame(checks)
print(checks.to_string(index=False))
checks.to_csv(RESULTS / "step0_centred_check.csv", index=False)

ds = xr.Dataset({k: ("time", v.to_numpy()) for k, v in series.items()}, coords={"time": GRID})
ds.attrs = {"title": "Centred and causal series, 15 min UTC grid (ANALYSIS_PLAN.md § 1)",
            "units": "speed m d-1; uplift m; accel m d-2; level m a.s.l.; tpos degC; rain mm h-1",
            "source": "Sugiyama et al. 2026, doi:10.17632/8wvtxg53ry.1; processing v1.0.2 doi:10.5281/zenodo.23257925"}
ds.to_netcdf(RESULTS / "series_15min.nc")

# %% [markdown]
# Peak shift between centred and causal processing, Period II (2–6 Jan).

# %%
P2 = ("2022-01-02", "2022-01-06")
BASE_VARS = ["tpos", "rain", "level", "uplift_GNSS1", "speed_GNSS1", "speed_GNSS2", "accel_GNSS1"]
rows = []
for var in BASE_VARS:
    for h in BANDWIDTHS:
        t = {proc: series[f"{var}_{proc}_h{h}"][P2[0]:P2[1]].idxmax() for proc in PROCS}
        rows.append({"series": var, "h": h, "peak_centred": t["centred"], "peak_causal": t["causal"],
                     "causal_minus_centred_h": (t["causal"] - t["centred"]) / pd.Timedelta("1h")})
peak_shift = pd.DataFrame(rows)
peak_shift.to_csv(RESULTS / "peak_shift.csv", index=False)
print(peak_shift.to_string(index=False))

# %% [markdown]
# ## Step 2 — onsets in Period II (and Period I, context only)

# %%
PERIODS = {
    "II": {"refs": {"R": ("2022-01-01 00:00", "2022-01-02 00:00"), "R'": ("2022-01-01 12:00", "2022-01-02 00:00")},
           "search": ("2022-01-02 00:00", "2022-01-06 00:00"), "vars": BASE_VARS},
    "I": {"refs": {"R": ("2021-12-20 00:00", "2021-12-21 00:00")},
          "search": ("2021-12-21 00:00", "2021-12-26 00:00"),
          "vars": ["tpos", "rain", "uplift_GNSS1", "speed_GNSS1", "speed_GNSS2"]},
}
GNSS1_GAP_I = (pd.Timestamp("2021-12-23 22:45"), pd.Timestamp("2021-12-25 12:20"))
RULES = ["T2", "T3", "T5", "CP"]


def onsets_for(s: pd.Series, refwin: tuple[str, str], search: tuple[str, str]) -> list[dict]:
    out = []
    for k in (2, 3, 5):
        out.append({"rule": f"T{k}", "onset": onset_threshold(s, refwin, search, k)})
    sv = s.dropna()
    cp = onset_changepoint(s, refwin[0], search)
    peak_t = s[search[0]:search[1]].idxmax() if s[search[0]:search[1]].notna().any() else None
    if cp and peak_t is not None and sv[refwin[0]:].index[0] > peak_t - pd.Timedelta("12h"):
        cp = None  # < 12 h of data before the maximum (GNSS2, Period I)
    out.append({"rule": "CP", "onset": cp[0] if cp else None, "lo": cp[1] if cp else None, "hi": cp[2] if cp else None})
    out.append({"rule": "peak", "onset": peak_t})
    return out


rows = []
for period, cfg in PERIODS.items():
    for var in cfg["vars"]:
        for proc in PROCS:
            for h in BANDWIDTHS:
                s = series[f"{var}_{proc}_h{h}"]
                for rname, refwin in cfg["refs"].items():
                    for r in onsets_for(s, refwin, cfg["search"]):
                        r.update(period=period, series=var, proc=proc, h=h, ref=rname)
                        o = r["onset"]
                        r["before_S"] = o is not None and o < pd.Timestamp(cfg["search"][0])
                        r["near_gap"] = (o is not None and var.endswith("GNSS1") and period == "I"
                                         and GNSS1_GAP_I[0] - pd.Timedelta("1h") <= o <= GNSS1_GAP_I[1] + pd.Timedelta("1h"))
                        rows.append(r)
onsets = pd.DataFrame(rows)[["period", "series", "proc", "h", "ref", "rule", "onset", "lo", "hi", "before_S", "near_gap"]]
raw_rain = aws.rain_intensity.to_series()["2022-01-02":"2022-01-06"]
raw_rain_onset = raw_rain[raw_rain > 0].index.min()
onsets = pd.concat([onsets, pd.DataFrame([{"period": "II", "series": "rain", "proc": "raw", "h": 0, "ref": "-",
                                           "rule": "first>0", "onset": raw_rain_onset}])], ignore_index=True)
onsets.to_csv(RESULTS / "onsets.csv", index=False)
print(f"raw rain onset (first rain_intensity > 0 after 2 Jan): {raw_rain_onset}")
print(onsets[(onsets.period == "II") & (onsets.proc == "causal") & (onsets.ref == "R")]
      .pivot_table(index=["series", "h"], columns="rule", values="onset", aggfunc="first").to_string())

# %% [markdown]
# Leads and the pre-registered labels (rules T2, T3, T5, CP × h × {R, R′}; peak reported
# separately and not part of the label).

# %%
RES_MIN = {"level": 0.5, "tpos": 5.0, "rain": 5.0}  # timing resolution, minutes; GNSS 7.5


def res_min(var: str) -> float:
    return RES_MIN.get(var, 7.5)


PAIRS = [(b, e) for b in ["level", "uplift_GNSS1"] for e in ["speed_GNSS1", "accel_GNSS1", "speed_GNSS2"]]
PAIRS += [(f, e) for f in ["tpos", "rain"] for e in ["speed_GNSS1", "accel_GNSS1", "speed_GNSS2", "level", "uplift_GNSS1"]]
on2 = onsets[onsets.period == "II"].set_index(["series", "proc", "h", "ref", "rule"]).onset
rows = []
for b, e in PAIRS:
    for proc in PROCS:
        for h in BANDWIDTHS:
            for rname in PERIODS["II"]["refs"]:
                for rule in RULES + ["peak"]:
                    ob, oe = on2.get((b, proc, h, rname, rule)), on2.get((e, proc, h, rname, rule))
                    lead = (oe - ob) / pd.Timedelta("1h") if pd.notna(ob) and pd.notna(oe) else np.nan
                    rows.append({"basal_or_forcing": b, "endpoint": e, "proc": proc, "h": h, "ref": rname,
                                 "rule": rule, "lead_h": lead,
                                 "resolution_h": np.hypot(res_min(b), res_min(e)) / 60})
leads = pd.DataFrame(rows)
leads.to_csv(RESULTS / "leads.csv", index=False)


def label(d: pd.DataFrame) -> pd.Series:
    n_all, v = len(d), d.lead_h.dropna()
    first = d.basal_or_forcing.iloc[0]
    if len(v) < n_all / 2:
        lab = "indeterminate"
    elif (v > 1).sum() > len(v) / 2:
        lab = f"{first} leads"
    elif (v.abs() <= 1).sum() > len(v) / 2:
        lab = "simultaneous"
    elif (v < -1).sum() > len(v) / 2:
        lab = f"{first} lags"
    else:
        lab = "indeterminate"
    return pd.Series({"n_combinations": n_all, "n_with_both_onsets": len(v),
                      "n_lead_gt_1h": int((v > 1).sum()), "n_abs_le_1h": int((v.abs() <= 1).sum()),
                      "n_lead_lt_-1h": int((v < -1).sum()), "median_lead_h": v.median(),
                      "min_lead_h": v.min(), "max_lead_h": v.max(), "label": lab + " (n = 1 event, one site)"})


labels = (leads[leads.rule.isin(RULES)].groupby(["basal_or_forcing", "endpoint", "proc"]).apply(label).reset_index())
labels.to_csv(RESULTS / "lead_labels.csv", index=False)
print(labels.to_string(index=False))

# %% [markdown]
# ## Step 3 — lead–lag structure (ℓ > 0: X leads Y)

# %%
WINDOWS = {"W1": ("2022-01-02 00:00", "2022-01-06 00:00"), "W2": ("2022-01-01 00:00", "2022-01-12 08:45"),
           "W3": ("2022-01-01 00:00", "2022-01-23 13:00"), "W4": ("2022-01-06 12:00", "2022-01-12 08:45")}
XC_PAIRS = [(x, y) for x in ["level", "uplift_GNSS1", "tpos", "rain"] for y in ["speed_GNSS1", "speed_GNSS2", "accel_GNSS1"]]
XC_PAIRS += [(x, y) for x in ["tpos", "rain"] for y in ["level", "uplift_GNSS1"]]
MAX_LAG = int(24 / STEP_H)
rows = []
for x, y in XC_PAIRS:
    for proc in PROCS:
        for h in BANDWIDTHS:
            for w, win in WINDOWS.items():
                X, Y = series[f"{x}_{proc}_h{h}"], series[f"{y}_{proc}_h{h}"]
                r = xcorr(X, Y, win, MAX_LAG)
                n = int((X[win[0]:win[1]].notna() & Y[win[0]:win[1]].notna()).sum())
                best = r.idxmax() if r.notna().any() else np.nan
                rows.append({"X": x, "Y": y, "proc": proc, "h": h, "window": w, "n_pairs_lag0": n,
                             "lag_star_h": best * STEP_H if pd.notna(best) else np.nan,
                             "r_star": r.max(), "r0": r.get(0)})
xc = pd.DataFrame(rows)
xc.to_csv(RESULTS / "xcorr.csv", index=False)
print(xc[(xc.proc == "causal") & xc.X.isin(["level", "uplift_GNSS1"])]
      .pivot_table(index=["X", "Y", "window"], columns="h", values="lag_star_h").to_string())

# %% [markdown]
# ## Step 4 — descriptive nested models (not held-out prediction across events)

# %%
LAGS_H = [0, 3, 6]
MODELS = {"M1": ["v"], "M2": ["v", "tpos", "rain"], "M3": ["v", "tpos", "rain", "level", "uplift_GNSS1"]}


def design(target_st: str, h: int, H: int) -> pd.DataFrame:
    hourly_idx = GRID[GRID.minute == 0]
    src = {"v": series[f"speed_{target_st}_causal_h{h}"]}
    for k in ["tpos", "rain", "level", "uplift_GNSS1"]:
        src[k] = series[f"{k}_causal_h{h}"]
    cols = {}
    for k, s in src.items():
        for lag in LAGS_H:
            cols[f"{k}_l{lag}"] = s.shift(int(lag / STEP_H)).reindex(hourly_idx)
    v = src["v"]
    cols["y"] = (v.shift(-int(H / STEP_H)) - v).reindex(hourly_idx)
    return pd.DataFrame(cols).dropna()


def xcols(m: str) -> list[str]:
    return [f"{k}_l{lag}" for k in MODELS[m] for lag in LAGS_H]


def fit(d: pd.DataFrame, m: str):
    return sm.OLS(d.y, sm.add_constant(d[xcols(m)], has_constant="add")).fit()


def predict(res, d: pd.DataFrame, m: str) -> np.ndarray:
    return res.predict(sm.add_constant(d[xcols(m)], has_constant="add"))


def rmse(e: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(e)))) if len(e) else np.nan


WIN_FIT = {"GNSS1": WINDOWS["W2"], "GNSS2": WINDOWS["W3"]}
TEST_4C = ("2022-01-01 00:00", "2022-01-06 12:00")
rows = []
for st in ["GNSS1", "GNSS2"]:
    for h in BANDWIDTHS:
        cp = onsets.set_index(["period", "series", "proc", "h", "ref", "rule"]).onset
        e_on = cp.get(("II", f"speed_{st}", "causal", h, "R", "CP"))
        e_pk = cp.get(("II", f"speed_{st}", "causal", h, "R", "peak"))
        for H in [1, 3, 6]:
            d = design(st, h, H)
            dw = d[WIN_FIT[st][0]:WIN_FIT[st][1]]
            # 4a: in-sample, same rows for all models
            for m in MODELS:
                r = fit(dw, m)
                rows.append({"part": "4a", "station": st, "h": h, "H": H, "model": m, "n_train": len(dw),
                             "r2": r.rsquared, "r2_adj": r.rsquared_adj, "aic": r.aic, "bic": r.bic})
            # 4b: rolling origin, expanding window
            preds = {m: [] for m in MODELS}
            truth, origins = [], []
            for i, t in enumerate(dw.index):
                train = dw[dw.index <= t - pd.Timedelta(hours=H)]
                if len(train) < 80:
                    continue
                origins.append(t)
                truth.append(dw.y.iloc[i])
                for m in MODELS:
                    preds[m].append(float(predict(fit(train, m), dw.iloc[[i]], m).iloc[0]))
            truth = np.array(truth)
            origins = pd.DatetimeIndex(origins)
            in_p2 = (origins >= P2[0]) & (origins < "2022-01-06")
            for m in MODELS:
                e = np.array(preds[m]) - truth
                rows.append({"part": "4b", "station": st, "h": h, "H": H, "model": m,
                             "first_origin": origins.min() if len(origins) else None, "n_test": len(e),
                             "event_onset_cp": e_on, "onset_before_first_origin": bool(len(origins) and e_on is not None and e_on < origins.min()),
                             "rmse": rmse(e), "rmse_period2": rmse(e[in_p2]), "n_test_period2": int(in_p2.sum())})
            # 4c: event hold-out, train after the event
            train = d["2022-01-06 12:00":WIN_FIT[st][1]]
            train = train[train.index <= pd.Timestamp(WIN_FIT[st][1]) - pd.Timedelta(hours=H)]
            test = d[TEST_4C[0]:TEST_4C[1]]
            onset_part = test[(e_on - pd.Timedelta("12h")):e_pk] if e_on is not None and e_pk is not None else test.iloc[:0]
            for m in MODELS:
                r = fit(train, m)
                rows.append({"part": "4c", "station": st, "h": h, "H": H, "model": m, "n_train": len(train),
                             "n_test": len(test), "rmse": rmse(predict(r, test, m) - test.y),
                             "n_test_onset": len(onset_part),
                             "rmse_onset": rmse(predict(r, onset_part, m) - onset_part.y) if len(onset_part) else np.nan})
nm = pd.DataFrame(rows)
for part in ["4b", "4c"]:
    for col in ["rmse", "rmse_period2", "rmse_onset"]:
        if col not in nm:
            continue
        sel = nm.part == part
        key = nm.loc[sel, ["station", "h", "H"]].astype(str).agg("|".join, axis=1)
        base = nm.loc[sel].assign(k=key).set_index(["k", "model"])[col]
        for ref_m in ["M1", "M2"]:
            nm.loc[sel, f"skill_{col}_vs_{ref_m}"] = [
                1 - base[(k, m)] / base[(k, ref_m)] for k, m in zip(key, nm.loc[sel, "model"])]
nm.to_csv(RESULTS / "nested_models.csv", index=False)

d4 = nm[(nm.part == "4c") & (nm.station == "GNSS1") & (nm.h == 3) & (nm.model == "M3")].set_index("H")
s_m3 = d4["skill_rmse_vs_M2"]
adds = s_m3.loc[3] > 0.10 and all(np.sign(s_m3.loc[[1, 6]]) == np.sign(s_m3.loc[3]))
print("4c skill of M3 vs M2 (GNSS1, h = 3 h) by horizon:", s_m3.round(3).to_dict())
print("Pre-registered label:", ("M3 adds information here" if adds else "M3 does not add information here"),
      "(n = 1 event, one site; not a held-out test across events)")
pd.DataFrame([{"rule": "M3 vs M2, 4c, GNSS1, h=3", **{f"skill_H{H}": s_m3.loc[H] for H in [1, 3, 6]},
               "label": ("M3 adds information here" if adds else "M3 does not add information here") + " (n = 1)"}]
             ).to_csv(RESULTS / "nested_models_label.csv", index=False)
print(nm[nm.part != "4a"].round(4).to_string(index=False))
