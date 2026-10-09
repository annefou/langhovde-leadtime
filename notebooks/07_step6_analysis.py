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
# # 07 — Step 6.2–6.5: onsets and leads from the kinematic GNSS (n = 1 event)
#
# `ANALYSIS_PLAN.md` § 6, pre-registered before processing. Inputs:
# - `results/kinematic_5min.nc` (notebook 06);
# - BH2201 water level with the measured on-glacier air-pressure correction (ADS
#   A20220506-001);
# - AWS forcing as in Step 1.
#
# **Still one event at one site.** Lead = onset(endpoint) − onset(basal), positive =
# basal first; the Step 6 label subtracts the method's own detection delay d(h).

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, ".")
from leadtime import (backward_rate, onset_changepoint, onset_threshold,  # noqa: E402
                      robust_causal_on_grid, smooth_on_grid, synthetic_onset_delay, xcorr)

ROOT = Path("..")
CLEAN, RESULTS = ROOT / "data" / "clean", ROOT / "results"
H6 = [0.25, 0.5, 1, 3]
PROCS = ["centred", "causal", "causal-robust"]
STEP_H = 0.25
GRID = pd.date_range("2021-12-31 00:00", "2022-01-07 00:00", freq="15min", inclusive="left")

# %% [markdown]
# ## 6.2 BH2201 level with measured air pressure

# %%
p = xr.open_dataset(CLEAN / "pressure_BH2201.nc").pressure.to_series().dropna()
aws_n = pd.read_csv(ROOT / "data/raw/nipr/A20220506-001/aws_220206.txt", parse_dates=["TIMESTAMP"],
                    index_col="TIMESTAMP", na_values=["NAN"])
pair = aws_n.AirPress.reindex(aws_n.AirPress.index.union(p.index)).interpolate("time").reindex(p.index)
pair_mean = pair.mean()
level_v102 = ((p - 9.2) * 1.0197 - 436.1)["2022-01-01":]
level_corr = ((p - 9.2 - (pair - pair_mean) / 100) * 1.0197 - 436.1)["2022-01-01":]
print(f"air pressure over BH2201 record: mean {pair_mean:.1f} hPa, range {pair.min():.1f}–{pair.max():.1f}")
print(f"Period II peak, corrected: {level_corr['2022-01-02':'2022-01-06'].max():.2f} m at "
      f"{level_corr['2022-01-02':'2022-01-06'].idxmax()}; uncorrected {level_v102['2022-01-02':'2022-01-06'].max():.2f} m")

# %% [markdown]
# ## 6.3 Series on the 15 min grid (kinematic GNSS; h = 0.25, 0.5, 1, 3 h)

# %%
kin = xr.open_dataset(RESULTS / "kinematic_5min.nc")
aws = xr.open_dataset(CLEAN / "aws.nc")
forcing = {"tpos": (aws.time.values, np.maximum(aws.air_temperature.values, 0.0)),
           "rain": (aws.time.values, aws.rain_intensity.values),
           "level": (level_corr.index.values, level_corr.to_numpy()),
           "level_v102": (level_v102.index.values, level_v102.to_numpy())}


def central_rate(x: np.ndarray, step_h: float) -> np.ndarray:
    out = np.full(x.size, np.nan)
    out[1:-1] = (x[2:] - x[:-2]) / (2 * step_h / 24)
    return out


series: dict[str, pd.Series] = {}
for h in H6:
    lag = int(max(1, h) / STEP_H)
    for st in ["GNSS1", "GNSS2"]:
        k = kin[[f"x_{st}_fixed", f"y_{st}_fixed", f"z_{st}_fixed"]].to_dataframe().dropna()
        t = k.index.values
        xyz = [k.iloc[:, i].to_numpy() - k.iloc[0, i] for i in range(3)]
        cen = [smooth_on_grid(t, v, GRID, h, causal=False) for v in xyz]
        cau = [smooth_on_grid(t, v, GRID, h, causal=True) for v in xyz]
        rob = [robust_causal_on_grid(t, v, GRID, h) for v in xyz]
        series[f"speed_{st}_centred_h{h}"] = pd.Series(np.hypot(central_rate(cen[0], STEP_H), central_rate(cen[1], STEP_H)), GRID)
        series[f"speed_{st}_causal_h{h}"] = pd.Series(np.hypot(backward_rate(cau[0], STEP_H), backward_rate(cau[1], STEP_H)), GRID)
        series[f"speed_{st}_causal-robust_h{h}"] = pd.Series(np.hypot(rob[0][1], rob[1][1]), GRID)
        for proc, z in (("centred", cen[2]), ("causal", cau[2]), ("causal-robust", rob[2][0])):
            series[f"uplift_{st}_{proc}_h{h}"] = pd.Series(z, GRID)
        for proc in PROCS:
            v = series[f"speed_{st}_{proc}_h{h}"].to_numpy()
            series[f"accel_{st}_{proc}_h{h}"] = pd.Series(backward_rate(v, STEP_H, lag), GRID)
    for name, (t, v) in forcing.items():
        for proc in ("centred", "causal"):
            series[f"{name}_{proc}_h{h}"] = pd.Series(smooth_on_grid(t, v, GRID, h, causal=proc == "causal"), GRID)
        series[f"{name}_causal-robust_h{h}"] = series[f"{name}_causal_h{h}"]
xr.Dataset({k: ("time", v.to_numpy()) for k, v in series.items()}, coords={"time": GRID}).to_netcdf(
    RESULTS / "step6_series_15min.nc")

# %% [markdown]
# ## 6.4 Detection delay at the kinematic noise level

# %%
noise = {}
for st in ["GNSS1", "GNSS2"]:
    y = kin[f"y_{st}_fixed"].sel(time=slice("2022-01-01", "2022-01-01 23:59")).to_series().dropna()
    tt = (y.index - y.index[0]) / pd.Timedelta("1D")
    noise[st] = float(np.std(y - np.polyval(np.polyfit(tt, y, 1), tt), ddof=1))
print("5 min kinematic northing noise SD (detrended, 1 Jan):", {k: round(v, 4) for k, v in noise.items()})
rows = []
for h in H6:
    for proc in ("causal", "causal-robust"):
        d = synthetic_onset_delay(h, noise["GNSS1"], 5, proc)
        rows.append({"h": h, "proc": proc, "noise_sd_m": noise["GNSS1"], "d_median_h": np.nanmedian(d),
                     "d_iqr_h": float(np.subtract(*np.nanpercentile(d, [75, 25]))), "n_defined": int(np.isfinite(d).sum())})
delay6 = pd.DataFrame(rows)
delay6.to_csv(RESULTS / "step6_delay_check.csv", index=False)
print(delay6.round(2).to_string(index=False))
D = {(r.proc, r.h): r.d_median_h for r in delay6.itertuples()} | {("centred", h): 0.0 for h in H6}

# %% [markdown]
# ## Onsets (same rules and windows as § 2)

# %%
REFS = {"R": ("2022-01-01 00:00", "2022-01-02 00:00"), "R'": ("2022-01-01 12:00", "2022-01-02 00:00")}
SEARCH = ("2022-01-02 00:00", "2022-01-06 00:00")
VARS = ["tpos", "rain", "level", "level_v102", "uplift_GNSS1", "speed_GNSS1", "speed_GNSS2", "accel_GNSS1"]
RULES = ["T2", "T3", "T5", "CP"]
rows = []
for var in VARS:
    for proc in PROCS:
        for h in H6:
            s = series[f"{var}_{proc}_h{h}"]
            for rname, ref in REFS.items():
                for k in (2, 3, 5):
                    rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": f"T{k}",
                                 "onset": onset_threshold(s, ref, SEARCH, k)})
                cp = onset_changepoint(s, ref[0], SEARCH)
                rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": "CP",
                             "onset": cp[0] if cp else None, "lo": cp[1] if cp else None, "hi": cp[2] if cp else None})
                w = s[SEARCH[0]:SEARCH[1]].dropna()
                rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": "peak",
                             "onset": w.idxmax() if len(w) else None})
on6 = pd.DataFrame(rows)
on6.to_csv(RESULTS / "step6_onsets.csv", index=False)
print(on6[(on6.ref == "R") & on6.rule.isin(["CP", "T3", "peak"])].pivot_table(
    index=["series", "rule"], columns=["proc", "h"], values="onset", aggfunc="first").to_string())

# %% [markdown]
# ## 6.5 Leads and the pre-registered, delay-corrected labels

# %%
idx = on6.set_index(["series", "proc", "h", "ref", "rule"]).onset
PAIRS = [(b, e) for b in ["level", "level_v102", "uplift_GNSS1"] for e in ["speed_GNSS1", "accel_GNSS1", "speed_GNSS2"]]
PAIRS += [(f, "level") for f in ["rain", "tpos"]]
rows = []
for b, e in PAIRS:
    for proc in PROCS:
        for h in H6:
            for rname in REFS:
                for rule in RULES + ["peak"]:
                    ob, oe = idx.get((b, proc, h, rname, rule)), idx.get((e, proc, h, rname, rule))
                    lead = (oe - ob) / pd.Timedelta("1h") if pd.notna(ob) and pd.notna(oe) else np.nan
                    d = D[(proc, h)] if e.startswith(("speed", "accel")) else 0.0
                    rows.append({"basal": b, "endpoint": e, "proc": proc, "h": h, "ref": rname, "rule": rule,
                                 "lead_h": lead, "d_h": d, "lead_minus_d_h": lead - d})
leads6 = pd.DataFrame(rows)
leads6.to_csv(RESULTS / "step6_leads.csv", index=False)


def label(v: pd.Series, n_all: int, first: str) -> str:
    v = v.dropna()
    if len(v) < n_all / 2:
        return "indeterminate"
    if (v > 1).sum() > len(v) / 2:
        return f"{first} leads"
    if (v.abs() <= 1).sum() > len(v) / 2:
        return "simultaneous within resolution"
    if (v < -1).sum() > len(v) / 2:
        return f"{first} lags"
    return "indeterminate"


rows = []
for (b, e, proc), g in leads6[leads6.rule.isin(RULES)].groupby(["basal", "endpoint", "proc"]):
    for col, kind in (("lead_minus_d_h", "delay-corrected (primary)"), ("lead_h", "uncorrected (d = 0)")):
        v = g[col]
        rows.append({"basal": b, "endpoint": e, "proc": proc, "version": kind, "n_combinations": len(g),
                     "n_defined": int(v.notna().sum()), "median_h": v.median(), "min_h": v.min(), "max_h": v.max(),
                     "label": label(v, len(g), b) + " (n = 1 event, one site)"
                     + (" — uses future data" if proc == "centred" else ""),
                     "status": "pre-registered (§ 6.5)"})
lab6 = pd.DataFrame(rows)
lab6.to_csv(RESULTS / "step6_lead_labels.csv", index=False)
print(lab6[lab6.basal.isin(["level", "rain", "tpos"])].round(2).to_string(index=False))

# %% [markdown]
# ## Cross-correlation (ℓ > 0: X leads Y), Period II window and 1–7 Jan

# %%
rows = []
for x in ["level", "uplift_GNSS1"]:
    for y in ["speed_GNSS1", "speed_GNSS2", "accel_GNSS1"]:
        for proc in PROCS:
            for h in H6:
                for w, win in {"W1": SEARCH, "W1-7": ("2022-01-01 00:00", "2022-01-06 23:45")}.items():
                    r = xcorr(series[f"{x}_{proc}_h{h}"], series[f"{y}_{proc}_h{h}"], win, 96)
                    best = r.idxmax() if r.notna().any() else np.nan
                    rows.append({"X": x, "Y": y, "proc": proc, "h": h, "window": w,
                                 "lag_star_h": best * STEP_H if pd.notna(best) else np.nan, "r_star": r.max(), "r0": r.get(0)})
xc6 = pd.DataFrame(rows)
xc6.to_csv(RESULTS / "step6_xcorr.csv", index=False)
print(xc6[(xc6.X == "level") & (xc6.Y == "speed_GNSS1")].pivot_table(
    index=["proc", "window"], columns="h", values=["lag_star_h", "r_star"]).round(2).to_string())
