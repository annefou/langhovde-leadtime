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
# # 08 — Step 7: is the ~12 h variability GNSS multipath or tide? (n = 1 event)
#
# `ANALYSIS_PLAN.md` § 7, pre-registered before computing. 7.1 sidereal filter with a
# pre-event template; 7.2 tidal regression fitted before the event and tested after it;
# 7.3 the § 6.3–6.5 onset analysis rerun on the cleaned series (via notebook 07).

# %%
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer

sys.path.insert(0, ".")
from leadtime import smooth_on_grid  # noqa: E402

ROOT = Path("..")
WORK, RESULTS, CLEAN = ROOT / "data/work/kinematic", ROOT / "results", ROOT / "data/clean"
ROVERS = {"GNSS1": "LG5", "GNSS2": "LG4"}
DAYS = pd.date_range("2021-12-31", "2022-01-06", freq="D")
TEMPLATE_DAYS = [pd.Timestamp("2021-12-31"), pd.Timestamp("2022-01-01")]
T_REPEAT = {"primary": 86154.0, "sensitivity": 86164.0}
TO_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32737", always_xy=True)
GRID15 = pd.date_range("2021-12-31", "2022-01-07", freq="15min", inclusive="left")


def read_pos(path: Path) -> pd.DataFrame:
    rows = [ln.split()[:6] for ln in path.read_text().splitlines() if ln and not ln.startswith("%")]
    df = pd.DataFrame(rows, columns=["date", "time", "lat", "lon", "h", "Q"])
    df.index = pd.to_datetime(df.date + " " + df.time, format="%Y/%m/%d %H:%M:%S.%f")
    return df[["lat", "lon", "h", "Q"]].astype(float)


def pos30(st: str) -> pd.DataFrame:
    d = pd.concat([read_pos(WORK / f"{ROVERS[st]}_{day:%Y%m%d}.pos") for day in DAYS]).sort_index()
    d = d[(d.Q == 1) & ~d.index.duplicated()]
    x, y = TO_UTM.transform(d.lon.to_numpy(), d.lat.to_numpy())
    d = pd.DataFrame({"x": x, "y": y, "z": d.h.to_numpy()}, index=d.index)
    g = d.groupby(d.index.floor("30s"))
    return g.median()[g.size() >= 15]


pos = {st: pos30(st) for st in ROVERS}
print({st: len(p) for st, p in pos.items()})

# %% [markdown]
# ## 7.1 Sidereal filter, template from 31 Dec and 1 Jan only

# %%
def residuals(p: pd.DataFrame, day: pd.Timestamp) -> pd.DataFrame:
    d = p[day: day + pd.Timedelta("1D") - pd.Timedelta("1s")]
    t = (d.index - day) / pd.Timedelta("1D")
    return d.apply(lambda c: c - np.polyval(np.polyfit(t, c, 1), t))


def sidereal_filter(p: pd.DataFrame, T: float) -> pd.DataFrame:
    tgt = p["2022-01-01":]
    ts = (tgt.index - tgt.index[0]).total_seconds().to_numpy()
    acc = np.zeros((len(tgt), 3))
    cnt = np.zeros(len(tgt))
    for tday in TEMPLATE_DAYS:
        r = residuals(p, tday)
        for n in range(1, 7):
            src_t = ((r.index + pd.Timedelta(seconds=n * T)) - tgt.index[0]).total_seconds().to_numpy()
            # only template days strictly before the target day
            ok_day = (tgt.index.floor("D") > tday).to_numpy()
            inside = (ts >= src_t.min()) & (ts <= src_t.max()) & ok_day
            if not inside.any():
                continue
            vals = np.column_stack([np.interp(ts[inside], src_t, r[c].to_numpy()) for c in "xyz"])
            acc[inside] += vals
            cnt[inside] += 1
    tmpl = np.where(cnt[:, None] > 0, acc / np.maximum(cnt, 1)[:, None], 0.0)
    out = p.copy()
    out.loc[tgt.index, ["x", "y", "z"]] = tgt[["x", "y", "z"]].to_numpy() - tmpl
    return out


def to5min(p: pd.DataFrame, st: str) -> dict[str, pd.Series]:
    g = p.groupby(p.index.floor("5min"))
    m = g.median()[g.size() >= 5]
    return {f"{c}_{st}_fixed": m[c] for c in "xyz"}


grid5 = pd.date_range("2021-12-31", "2022-01-07", freq="5min", inclusive="left")
filtered = {}
for tag, T in T_REPEAT.items():
    data = {}
    for st in ROVERS:
        f = sidereal_filter(pos[st], T)
        filtered[tag, st] = f
        data |= to5min(f, st)
    ds = xr.Dataset({k: ("time", v.reindex(grid5).to_numpy()) for k, v in data.items()}, coords={"time": grid5})
    ds.attrs = {"title": f"Kinematic GNSS 5 min, sidereal-filtered (T = {T:.0f} s, template 31 Dec + 1 Jan)"}
    ds.to_netcdf(RESULTS / f"kinematic_5min_sidereal_{tag}.nc")
raw5 = xr.open_dataset(RESULTS / "kinematic_5min.nc")

# %% [markdown]
# Band variance (10–15 h) of centred h = 1 h speed, 2–6 Jan, before vs after.

# %%
def speed_h1(ds: xr.Dataset, st: str) -> pd.Series:
    k = ds[[f"x_{st}_fixed", f"y_{st}_fixed"]].to_dataframe().dropna()
    t = k.index.values
    sm = [smooth_on_grid(t, k.iloc[:, i].to_numpy() - k.iloc[0, i], GRID15, 1.0, causal=False) for i in range(2)]
    rate = [np.r_[np.nan, (v[2:] - v[:-2]) / (0.5 / 24), np.nan] for v in sm]
    return pd.Series(np.hypot(*rate), GRID15)


def band_var(s: pd.Series, lo_h: float = 10, hi_h: float = 15) -> float:
    v = s["2022-01-02":"2022-01-06 23:45"].interpolate().bfill().ffill().to_numpy()
    v = v - v.mean()
    f = np.fft.rfftfreq(v.size, 0.25)
    pw = np.abs(np.fft.rfft(v)) ** 2
    sel = (f >= 1 / hi_h) & (f <= 1 / lo_h)
    return float(pw[sel].sum() / v.size ** 2 * 2)


rows = []
spd = {}
for st in ROVERS:
    spd["raw", st] = speed_h1(raw5, st)
    for tag in T_REPEAT:
        spd[tag, st] = speed_h1(xr.open_dataset(RESULTS / f"kinematic_5min_sidereal_{tag}.nc"), st)
    b0 = band_var(spd["raw", st])
    for tag in T_REPEAT:
        b1 = band_var(spd[tag, st])
        rows.append({"station": st, "T_repeat": tag, "band_var_raw": b0, "band_var_filtered": b1,
                     "reduction": 1 - b1 / b0,
                     "sd_speed_raw": float(spd["raw", st]["2022-01-02":"2022-01-06"].std()),
                     "sd_speed_filtered": float(spd[tag, st]["2022-01-02":"2022-01-06"].std())})
sid = pd.DataFrame(rows)
red = sid[sid.T_repeat == "primary"].set_index("station").reduction
label71 = ("multipath dominant" if (red >= 0.5).all() else "multipath minor" if (red < 0.2).all() else "mixed")
sid["label"] = label71 + " (pre-registered § 7.1)"
sid.to_csv(RESULTS / "step7_sidereal.csv", index=False)
print(sid.round(4).to_string(index=False))

# %% [markdown]
# ## 7.2 Tidal regression: fit 31 Dec 12:00 – 2 Jan 00:00, held-out R² on 5–6 Jan

# %%
tide = xr.open_dataset(CLEAN / "tide.nc").tide.to_series()
tide.index = tide.index - pd.Timedelta("2.75h")
eta = tide.resample("15min").mean()
eta = eta - eta.mean()
deta = eta.diff() / (0.25 / 24)
FIT, TEST = ("2021-12-31 12:00", "2022-01-01 23:45"), ("2022-01-05 00:00", "2022-01-06 23:45")
rows, model = [], {}
for st in ROVERS:
    v = spd["primary", st]
    best = None
    for lag_steps in range(-24, 25):
        X = pd.DataFrame({"eta": eta.shift(lag_steps), "deta": deta.shift(lag_steps)}).reindex(GRID15)
        d = pd.concat([v.rename("v"), X], axis=1).dropna()
        tr, te = d[FIT[0]:FIT[1]], d[TEST[0]:TEST[1]]
        A = np.column_stack([np.ones(len(tr)), tr.eta, tr.deta])
        coef, *_ = np.linalg.lstsq(A, tr.v, rcond=None)
        r2_fit = 1 - ((tr.v - A @ coef) ** 2).sum() / ((tr.v - tr.v.mean()) ** 2).sum()
        if best is None or r2_fit > best["r2_fit"]:
            pred = coef[0] + coef[1] * te.eta + coef[2] * te.deta
            r2_te = 1 - ((te.v - pred) ** 2).sum() / ((te.v - te.v.mean()) ** 2).sum()
            best = {"station": st, "lag_h": lag_steps * 0.25, "a": coef[0], "b": coef[1], "c": coef[2],
                    "r2_fit": r2_fit, "r2_heldout": r2_te, "n_fit": len(tr), "n_test": len(te)}
    rows.append(best)
    model[st] = {"lag_h": best["lag_h"], "b": best["b"], "c": best["c"]}
tid = pd.DataFrame(rows)
r2 = tid.set_index("station").r2_heldout
label72 = ("tide explains the sub-daily variability" if (r2 >= 0.25).all() else "weak" if (r2 < 0.10).all() else "mixed")
tid["label"] = label72 + " (pre-registered § 7.2)"
tid.to_csv(RESULTS / "step7_tide.csv", index=False)
(RESULTS / "step7_tide_model.json").write_text(json.dumps(model, indent=1))
print(tid.round(4).to_string(index=False))

# %% [markdown]
# ## 7.3 Onset analysis on (a) sidereal-filtered and (b) sidereal + tide-corrected series

# %%
for prefix, extra in (("step7a", {}), ("step7b", {"STEP6_TIDE_MODEL": str((RESULTS / "step7_tide_model.json").resolve())})):
    env = os.environ | {"STEP6_KIN": "kinematic_5min_sidereal_primary.nc", "STEP6_PREFIX": prefix,
                        "STEP6_NOISE_FROM_DIFF": "1"} | extra
    r = subprocess.run([sys.executable, "-W", "ignore", "07_step6_analysis.py"], env=env, capture_output=True, text=True)
    (RESULTS / "logs" / f"07_{prefix}.log").write_text(r.stdout + r.stderr)
    if r.returncode:
        raise RuntimeError(f"{prefix} failed; see results/logs/07_{prefix}.log")
    print(prefix, "done")
