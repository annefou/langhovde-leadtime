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
# # 09 — Step 8: common-mode test with rovers LG1 and LGLK (n = 1 event)
#
# `ANALYSIS_PLAN.md` § 8, pre-registered. LG1 and LGLK are processed exactly as in
# § 6.1 / A4: kinematic against LGFX, the same single base position and base ephemeris,
# fixed only, 5 min medians (raw format detected per file: LG1 is a Septentrio
# receiver). Their mean residual is the common mode removed from
# GNSS1 and GNSS2.

# %%
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

ROOT = Path("..").resolve()
BIN = ROOT / ".tools/rtklib/bin"
RAW = ROOT / "data/raw/nipr/A20220506-004"
WORK, RESULTS = ROOT / "data/work/kinematic", ROOT / "results"
TO_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32737", always_xy=True)
NEW = {"LG1": ("LG1", "GPS1", pd.date_range("2021-12-31", "2022-01-06")),
       "LGLK": ("LGLK", "LGLK", pd.date_range("2021-12-31", "2022-01-03"))}
GRID5 = pd.date_range("2021-12-31", "2022-01-07", freq="5min", inclusive="left")
GRID15 = pd.date_range("2021-12-31", "2022-01-07", freq="15min", inclusive="left")


def raw_format(path: Path) -> str:
    """'sbf' for Septentrio (starts with '$@'), 'javad' for JAVAD JPS ('MF' or '~~')."""
    with open(path, "rb") as fh:
        head = fh.read(2)
    if head == b"$@":
        return "sbf"
    if head in (b"MF", b"~~"):  # JAVAD file header, or a split file starting with a JAVAD message
        return "javad"
    raise ValueError(f"unknown raw format: {path}")


def rinex(path: Path, fmt: str) -> tuple[Path, Path]:
    o, n = WORK / f"{path.name}.obs", WORK / f"{path.name}.nav"
    if not o.exists():
        subprocess.run([str(BIN / "convbin"), "-r", fmt, "-o", str(o), "-n", str(n), str(path)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return o, n


def base_ecef() -> list[str]:
    """Same base position as notebook 06: mean of LGFX single-point solutions, 31 Dec – 6 Jan."""
    xyz = []
    for f in sorted(WORK.glob("LGFX_*_single.pos")):
        if f.name.endswith("_events.pos"):
            continue
        xyz += [[float(v) for v in ln.split()[2:5]] for ln in f.read_text().splitlines() if ln and not ln.startswith("%")]
    return [str(v) for v in np.array(xyz).mean(axis=0)]


BASE = base_ecef()
print("base ECEF:", [round(float(v), 3) for v in BASE])


def process(st: str, day: pd.Timestamp) -> pd.DataFrame:
    folder, stem, _ = NEW[st]
    out = WORK / f"{st}_{day:%Y%m%d}.pos"
    if not out.exists():
        raws = sorted((RAW / folder).glob(f"{stem}{day.dayofyear:03d}?.{day.year % 100:02d}"))
        b_obs, b_nav = rinex(RAW / "LGFX" / f"LGFX{day.dayofyear:03d}0.{day.year % 100:02d}", "sbf")
        text = ""
        for i, r in enumerate(raws):
            ro, _ = rinex(r, raw_format(r))
            part = WORK / f"{st}_{day:%Y%m%d}_{i}.pos"
            subprocess.run([str(BIN / "rnx2rtkp"), "-p", "2", "-m", "15", "-sys", "G", "-f", "2", "-v", "3",
                            "-u", "-t", "-d", "1", "-r", *BASE, "-ts", f"{day:%Y/%m/%d}", "00:00:00",
                            "-te", f"{day:%Y/%m/%d}", "23:59:59", "-o", str(part), str(ro), str(b_obs), str(b_nav)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            text += part.read_text()
        out.write_text(text)
    rows = [ln.split()[:6] for ln in out.read_text().splitlines() if ln and not ln.startswith("%")]
    d = pd.DataFrame(rows, columns=["date", "time", "lat", "lon", "h", "Q"])
    d.index = pd.to_datetime(d.date + " " + d.time, format="%Y/%m/%d %H:%M:%S.%f")
    return d[["lat", "lon", "h", "Q"]].astype(float)


# %% [markdown]
# ## 8.1 Kinematic LG1 and LGLK; fix rates

# %%
pos5, fix_rows = {}, []
for st, (_, _, days) in NEW.items():
    parts = []
    for day in days:
        d = process(st, day)
        d = d[~d.index.duplicated()]
        rate = float((d.Q == 1).mean()) if len(d) else 0.0
        fix_rows.append({"station": st, "day": day.date(), "n_epochs": len(d), "fix_rate": rate})
        if rate >= 0.5:
            parts.append(d[d.Q == 1])
    d = pd.concat(parts).sort_index()
    x, y = TO_UTM.transform(d.lon.to_numpy(), d.lat.to_numpy())
    e = pd.DataFrame({"x": x, "y": y, "z": d.h.to_numpy()}, index=d.index)
    g = e.groupby(e.index.floor("5min"))
    pos5[st] = g.median()[g.size() >= 150].reindex(GRID5)
fixdf = pd.DataFrame(fix_rows)
fixdf.to_csv(RESULTS / "step8_fix_rate.csv", index=False)
print(fixdf.round(3).to_string(index=False))

kin = xr.open_dataset(RESULTS / "kinematic_5min.nc")
for st in ["GNSS1", "GNSS2"]:
    pos5[st] = pd.DataFrame({c: kin[f"{c}_{st}_fixed"].to_series() for c in "xyz"}).reindex(GRID5)

# %% [markdown]
# ## 8.2 Per-day detrended residuals; common mode = mean of LG1 and LGLK

# %%
def residuals(p: pd.DataFrame) -> pd.DataFrame:
    out = p.copy() * np.nan
    for day, d in p.groupby(p.index.floor("D")):
        ok = d.notna().all(axis=1)
        if ok.sum() < 50:
            continue
        t = ((d.index - day) / pd.Timedelta("1D")).to_numpy()
        for c in "xyz":
            coef = np.polyfit(t[ok], d[c][ok], 1)
            out.loc[d.index, c] = d[c] - np.polyval(coef, t)
    return out


def bandpass(s: pd.Series, lo_h: float = 6, hi_h: float = 18) -> pd.Series:
    v = s.interpolate(limit=12).to_numpy()
    ok = np.isfinite(v)
    w = np.where(ok, v, 0.0)
    f = np.fft.rfftfreq(w.size, 5 / 60)
    F = np.fft.rfft(w - w[ok].mean() * ok)
    F[(f < 1 / hi_h) | (f > 1 / lo_h)] = 0
    return pd.Series(np.where(ok, np.fft.irfft(F, w.size), np.nan), s.index)


res = {st: residuals(p) for st, p in pos5.items()}
cm = pd.concat([res["LG1"], res["LGLK"]]).groupby(level=0).mean().reindex(GRID5)

rows = []
names = ["GNSS1", "GNSS2", "LG1", "LGLK"]
for i, a in enumerate(names):
    for b in names[i + 1:]:
        for c, lab in zip("xyz", "ENU"):
            ba, bb = bandpass(res[a][c]), bandpass(res[b][c])
            ok = ba.notna() & bb.notna()
            rows.append({"pair": f"{a}-{b}", "component": lab, "n": int(ok.sum()),
                         "r_bandpassed_6_18h": float(np.corrcoef(ba[ok], bb[ok])[0, 1]) if ok.sum() > 50 else np.nan})
corr = pd.DataFrame(rows)
corr.to_csv(RESULTS / "step8_pair_correlation.csv", index=False)
print(corr.pivot(index="pair", columns="component", values="r_bandpassed_6_18h").round(2).to_string())

# %% [markdown]
# ## 8.3 Band variance of h = 1 h speed before vs after removing the common mode

# %%
def speed_h1(p: pd.DataFrame) -> pd.Series:
    d = p.dropna()
    t = d.index.values
    sm = [smooth_on_grid(t, d[c].to_numpy() - d[c].iloc[0], GRID15, 1.0, causal=False) for c in "xy"]
    rate = [np.r_[np.nan, (v[2:] - v[:-2]) / (0.5 / 24), np.nan] for v in sm]
    return pd.Series(np.hypot(*rate), GRID15)


def band_var(s: pd.Series, lo_h: float = 10, hi_h: float = 15) -> float:
    v = s["2022-01-02":"2022-01-06 23:45"].interpolate().bfill().ffill().to_numpy()
    v = v - v.mean()
    f = np.fft.rfftfreq(v.size, 0.25)
    pw = np.abs(np.fft.rfft(v)) ** 2
    return float(pw[(f >= 1 / hi_h) & (f <= 1 / lo_h)].sum() / v.size ** 2 * 2)


corrected, rows = {}, []
for st in ["GNSS1", "GNSS2"]:
    c = pos5[st] - cm.fillna(0.0)
    corrected[st] = c
    b0, b1 = band_var(speed_h1(pos5[st])), band_var(speed_h1(c))
    rows.append({"station": st, "band_var_raw": b0, "band_var_cm_removed": b1, "reduction": 1 - b1 / b0})
cmt = pd.DataFrame(rows)
red = cmt.set_index("station").reduction
label = ("common-mode dominant" if (red >= 0.5).all() else "common-mode minor" if (red < 0.2).all() else "mixed")
cmt["label"] = label + " (pre-registered § 8.3)"
cmt.to_csv(RESULTS / "step8_common_mode.csv", index=False)
print(cmt.round(4).to_string(index=False))

ds = xr.Dataset({f"{c}_{st}_fixed": ("time", corrected[st][c].to_numpy()) for st in corrected for c in "xyz"},
                coords={"time": GRID5})
ds.attrs = {"title": "Kinematic GNSS 5 min, common mode (mean of LG1, LGLK residuals) removed"}
ds.to_netcdf(RESULTS / "kinematic_5min_cm.nc")

# %% [markdown]
# ## 8.4 Onset analysis on the corrected positions (unless common mode is minor)

# %%
if label != "common-mode minor":
    env = os.environ | {"STEP6_KIN": "kinematic_5min_cm.nc", "STEP6_PREFIX": "step8", "STEP6_NOISE_FROM_DIFF": "1"}
    r = subprocess.run([sys.executable, "-W", "ignore", "07_step6_analysis.py"], env=env, capture_output=True, text=True)
    (RESULTS / "logs" / "07_step8.log").write_text(r.stdout + r.stderr)
    if r.returncode:
        raise RuntimeError("step8 onset analysis failed; see results/logs/07_step8.log")
    print("step8 onset analysis done")
else:
    print("common mode minor: § 6.3–6.5 not rerun (pre-registered § 8.4)")
