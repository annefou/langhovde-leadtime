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
# # 10 — Step 9: bed tilt in borehole BH2201 (n = 1 event)
#
# `ANALYSIS_PLAN.md` § 9, pre-registered after a format-only inspection. Data: ADS
# A20220506-002 (field "borehole 3" = paper BH2201), Sugiyama, Minowa, Kondo & Aoki
# (2022), CC BY 4.0. 100 Hz 3-axis accelerometer → 1 min medians → tilt (µrad).
# Pressure and GNSS1 speed come from the Step 6 series (`results/step6_series_15min.nc`).

# %%
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, ".")
from leadtime import backward_rate, onset_changepoint, onset_threshold, smooth_on_grid  # noqa: E402

ROOT = Path("..")
RAW = ROOT / "data/raw/nipr/A20220506-002"
RESULTS = ROOT / "results"
DAYS = ["01211231"] + [f"012201{d:02d}" for d in range(1, 7)]
H = [0.25, 0.5, 1, 3]
PROCS = ["centred", "causal"]
STEP_H = 0.25
GRID = pd.date_range("2021-12-31 00:00", "2022-01-07 00:00", freq="15min", inclusive="left")
REFS = {"R": ("2022-01-01 00:00", "2022-01-02 00:00"), "R'": ("2022-01-01 12:00", "2022-01-02 00:00")}
SEARCH = ("2022-01-02 00:00", "2022-01-06 00:00")

# %% [markdown]
# ## 9.1 Accelerometer → 1 min medians; daily GPS clock re-sync steps

# %%
def read10(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, header=None, names=["t", "a1", "a2", "a3", "c4"], skipinitialspace=True)
    day = pd.Timestamp(path.stem[:8])
    t = pd.to_timedelta(d.t)
    d.index = day + t
    return d[["a1", "a2", "a3", "c4"]]


parts = []
for day in DAYS:
    for f in sorted((RAW / day).glob("*/2022*.txt")) + sorted((RAW / day).glob("*/2021*.txt")):
        d = read10(f)
        parts.append(d.groupby(d.index.floor("1min")).median())
acc = pd.concat(parts).sort_index()
acc = acc[~acc.index.duplicated()]
print(f"1 min accelerometer: {len(acc)} minutes, {acc.index[0]} – {acc.index[-1]}")

steps = []
for day in DAYS:
    log = RAW / day / "LOGFILE.TXT"
    if log.exists():
        for ln in log.read_text(errors="replace").splitlines():
            m = re.match(r"(\S+ \S+)\s+GPS Set Time OK\. Set Time : (\S+ \S+)\.", ln)
            if m:
                a, b = pd.Timestamp(m.group(1)), pd.Timestamp(m.group(2))
                steps.append({"logged": a, "set_to": b, "step_s": (b - a).total_seconds()})
steps = pd.DataFrame(steps)
steps.to_csv(RESULTS / "step9_clock_steps.csv", index=False)
print(steps.to_string(index=False))

# %%
g = np.sqrt(acc.a1 ** 2 + acc.a2 ** 2 + acc.a3 ** 2)
tilt = pd.DataFrame({"th1": np.arcsin(acc.a2 / g) * 1e6, "th2": np.arcsin(acc.a3 / g) * 1e6}, index=acc.index)
ref = tilt[REFS["R"][0]:REFS["R"][1]].mean()
tilt["change"] = np.hypot(tilt.th1 - ref.th1, tilt.th2 - ref.th2)
r = tilt[REFS["R"][0]:REFS["R"][1]]
tt = (r.index - r.index[0]) / pd.Timedelta("1h")
drift = {c: float(np.polyfit(tt, r[c], 1)[0]) for c in ("th1", "th2")}
print("tilt drift in R (µrad/h):", {k: round(v, 3) for k, v in drift.items()})
z = (tilt.change - tilt.change.rolling(61, center=True, min_periods=10).median()).abs()
spikes = int((z > 10 * z.std()).sum())
print("single-minute spikes > 10 SD:", spikes)
xr.Dataset({c: ("time", tilt[c].to_numpy()) for c in tilt}, coords={"time": tilt.index.values}).to_netcdf(
    RESULTS / "step9_tilt_1min.nc")

# %% [markdown]
# ## Series on the 15 min grid

# %%
s6 = xr.open_dataset(RESULTS / "step6_series_15min.nc")
series = {}
t = tilt.index.values
for h in H:
    for proc in PROCS:
        ch = smooth_on_grid(t, tilt.change.to_numpy(), GRID, h, causal=proc == "causal", max_gap_h=1)
        series[f"tilt_{proc}_h{h}"] = pd.Series(ch, GRID)
        series[f"tiltrate_{proc}_h{h}"] = pd.Series(backward_rate(ch, STEP_H, int(max(1, h) / STEP_H)), GRID)
        for v in ("level", "speed_GNSS1"):
            series[f"{v}_{proc}_h{h}"] = s6[f"{v}_{proc}_h{h}"].to_series().reindex(GRID)
xr.Dataset({k: ("time", v.to_numpy()) for k, v in series.items()}, coords={"time": GRID}).to_netcdf(
    RESULTS / "step9_series_15min.nc")

# %% [markdown]
# ## 9.2 Quiet-period check

# %%
def bandpass(s: pd.Series, lo_h: float = 6, hi_h: float = 18) -> pd.Series:
    v = s.interpolate().bfill().ffill().to_numpy()
    f = np.fft.rfftfreq(v.size, STEP_H)
    F = np.fft.rfft(v - v.mean())
    F[(f < 1 / hi_h) | (f > 1 / lo_h)] = 0
    return pd.Series(np.fft.irfft(F, v.size), s.index)


s1 = series["tilt_centred_h1"]
quiet_sd = float(bandpass(s1)["2022-01-01":"2022-01-01 23:45"].std())
amp = float(s1[SEARCH[0]:SEARCH[1]].max() - s1[SEARCH[0]:SEARCH[1]].min())
label92 = "tilt resolves the event" if amp >= 5 * quiet_sd else "tilt does not resolve the event"
q = pd.DataFrame([{"quiet_bandpassed_sd_urad": quiet_sd, "event_amplitude_urad": amp, "ratio": amp / quiet_sd,
                   "drift_R_th1_urad_per_h": drift["th1"], "drift_R_th2_urad_per_h": drift["th2"],
                   "spikes_gt_10sd": spikes, "label": label92 + " (pre-registered § 9.2)"}])
q.to_csv(RESULTS / "step9_quiet_check.csv", index=False)
print(q.round(3).to_string(index=False))

# %% [markdown]
# ## 9.3 Onsets and leads (§ 2 rules, § 6.5 decision rule)

# %%
delay = pd.read_csv(RESULTS / "step6_delay_check.csv")
D = {(r.proc, r.h): r.d_median_h for r in delay.itertuples()} | {("centred", h): 0.0 for h in H}
rows = []
for var in ("level", "tilt", "tiltrate", "speed_GNSS1"):
    for proc in PROCS:
        for h in H:
            s = series[f"{var}_{proc}_h{h}"]
            for rname, rw in REFS.items():
                for k in (2, 3, 5):
                    rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": f"T{k}",
                                 "onset": onset_threshold(s, rw, SEARCH, k)})
                cp = onset_changepoint(s, rw[0], SEARCH)
                rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": "CP", "onset": cp[0] if cp else None})
                w = s[SEARCH[0]:SEARCH[1]].dropna()
                rows.append({"series": var, "proc": proc, "h": h, "ref": rname, "rule": "peak",
                             "onset": w.idxmax() if len(w) else None})
on9 = pd.DataFrame(rows)
on9.to_csv(RESULTS / "step9_onsets.csv", index=False)
print(on9[(on9.ref == "R") & on9.rule.isin(["CP", "T3", "peak"])].pivot_table(
    index=["series", "rule"], columns=["proc", "h"], values="onset", aggfunc="first").to_string())

idx = on9.set_index(["series", "proc", "h", "ref", "rule"]).onset
rows = []
for b, e in [("level", "tilt"), ("level", "tiltrate"), ("tilt", "speed_GNSS1"), ("level", "speed_GNSS1")]:
    for proc in PROCS:
        for h in H:
            for rname in REFS:
                for rule in ["T2", "T3", "T5", "CP", "peak"]:
                    ob, oe = idx.get((b, proc, h, rname, rule)), idx.get((e, proc, h, rname, rule))
                    lead = (oe - ob) / pd.Timedelta("1h") if pd.notna(ob) and pd.notna(oe) else np.nan
                    d = D[(proc, h)] if e == "speed_GNSS1" else 0.0
                    rows.append({"first": b, "second": e, "proc": proc, "h": h, "ref": rname, "rule": rule,
                                 "lead_h": lead, "d_h": d, "lead_minus_d_h": lead - d})
leads9 = pd.DataFrame(rows)
leads9.to_csv(RESULTS / "step9_leads.csv", index=False)


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
for (b, e, proc), gq in leads9[leads9.rule != "peak"].groupby(["first", "second", "proc"]):
    v = gq.lead_minus_d_h
    rows.append({"first": b, "second": e, "proc": proc, "n_combinations": len(gq), "n_defined": int(v.notna().sum()),
                 "median_h": v.median(), "min_h": v.min(), "max_h": v.max(),
                 "label": label(v, len(gq), b) + " (n = 1 event, one site)" + (" — uses future data" if proc == "centred" else ""),
                 "status": "pre-registered (§ 9.3)" + ("" if label92.startswith("tilt resolves") else "; no label per § 9.2")})
lab9 = pd.DataFrame(rows)
lab9.to_csv(RESULTS / "step9_lead_labels.csv", index=False)
print(lab9.round(2).to_string(index=False))
