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
# # 06 — Kinematic GNSS from the raw 1 s data (Step 6.1)
#
# `ANALYSIS_PLAN.md` § 6.1, pre-registered before processing. Data: ADS A20220506-004,
# Sugiyama, Minowa, Kondo & Aoki (2022), CC BY 4.0. Rovers LG5 (GNSS1, grounded) and
# LG4 (GNSS2, afloat) against base LGFX (rock), 31 Dec 2021 – 6 Jan 2022.
#
# RTKLIB v2.5.1 (`scripts/build_rtklib.sh`):
# - `convbin -r javad` to RINEX;
# - `rnx2rtkp`: relative kinematic, GPS, L1 + L2, mask 15°, ratio 3, continuous
#   ambiguity resolution, base position = average of its single-point solutions;
#   output in UTC.
#
# The filter runs **forward only** (no `-c` combined solution), so each position uses
# only past and present epochs. It is processed day by day, so the filter restarts at
# 00:00 UTC.
#
# Output: `results/kinematic_5min.nc` (5 min medians of fixed 1 s epochs, UTM 37S) and
# the consistency check against the deposit's 15 min positions (stop rule).

# %%
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer

ROOT = Path("..").resolve()
BIN = ROOT / ".tools" / "rtklib" / "bin"
RAW = ROOT / "data" / "raw" / "nipr" / "A20220506-004"
WORK = ROOT / "data" / "work" / "kinematic"
RESULTS = ROOT / "results"
WORK.mkdir(parents=True, exist_ok=True)
if not (BIN / "rnx2rtkp").exists():
    subprocess.run(["bash", str(ROOT / "scripts" / "build_rtklib.sh")], check=True)

DAYS = pd.date_range("2021-12-31", "2022-01-06", freq="D")
STATIONS = {"LG5": ("LG5", "GPS5"), "LG4": ("LG4", "LG04"), "LGFX": ("LGFX", "LGFX")}
ROVERS = {"GNSS1": "LG5", "GNSS2": "LG4"}
TO_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32737", always_xy=True)


def raw_files(station: str, day: pd.Timestamp) -> list[Path]:
    folder, stem = STATIONS[station]
    pat = f"{stem}{day.dayofyear:03d}?.{day.year % 100:02d}"
    return sorted((RAW / folder).glob(pat))


def to_rinex(station: str, day: pd.Timestamp) -> tuple[list[Path], list[Path]]:
    obs, nav = [], []
    for f in raw_files(station, day):
        o, n = WORK / f"{f.name}.obs", WORK / f"{f.name}.nav"
        if not o.exists():
            subprocess.run([str(BIN / "convbin"), "-r", "javad", "-o", str(o), "-n", str(n), str(f)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        obs.append(o)
        nav.append(n)
    return obs, nav


# %% [markdown]
# ## Conversion and kinematic processing, one rover-day at a time

# %%
def process_day(rover: str, day: pd.Timestamp) -> Path | None:
    out = WORK / f"{rover}_{day:%Y%m%d}.pos"
    if out.exists():
        return out
    r_obs, r_nav = to_rinex(rover, day)
    b_obs, _ = to_rinex("LGFX", day)
    if not r_obs or not b_obs:
        print(f"{rover} {day:%Y-%m-%d}: missing rover or base files")
        return None
    # Split rover days (two raw files) are processed per file; each is its own run.
    parts = []
    for i, ro in enumerate(r_obs):
        part = WORK / f"{rover}_{day:%Y%m%d}_{i}.pos"
        cmd = [str(BIN / "rnx2rtkp"), "-p", "2", "-m", "15", "-sys", "G", "-f", "2", "-v", "3",
               "-u", "-t", "-d", "1",
               "-ts", f"{day:%Y/%m/%d}", "00:00:00",
               "-te", f"{day:%Y/%m/%d}", "23:59:59",
               "-o", str(part), str(ro), *map(str, b_obs), str(r_nav[i])]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        parts.append(part)
    out.write_text("".join(p.read_text() for p in parts))
    return out


def read_pos(path: Path) -> pd.DataFrame:
    rows = [ln.split() for ln in path.read_text().splitlines() if ln and not ln.startswith("%")]
    df = pd.DataFrame([r[:6] for r in rows], columns=["date", "time", "lat", "lon", "h", "Q"])
    df.index = pd.to_datetime(df.date + " " + df.time, format="%Y/%m/%d %H:%M:%S.%f")
    return df[["lat", "lon", "h", "Q"]].astype({"lat": float, "lon": float, "h": float, "Q": int})


sols = {}
for name, st in ROVERS.items():
    parts = [read_pos(p) for d in DAYS if (p := process_day(st, d)) is not None]
    sols[name] = pd.concat(parts).sort_index()
    sols[name] = sols[name][~sols[name].index.duplicated()]

# %% [markdown]
# ## Fix rate per day (Q = 1 fixed, Q = 2 float)

# %%
fix = pd.concat({n: s.Q.groupby(s.index.floor("D")).agg(
    n_epochs="size", fix_rate=lambda q: (q == 1).mean(), float_rate=lambda q: (q == 2).mean())
    for n, s in sols.items()}, names=["station", "day"])
fix.to_csv(RESULTS / "kinematic_fix_rate.csv")
print(fix.round(3).to_string())
low = fix.loc[(slice(None), slice("2022-01-02", "2022-01-04")), "fix_rate"] < 0.5
print("days 2–4 Jan with fix rate < 50 %:", list(low[low].index))

# %% [markdown]
# ## 5 min medians of fixed epochs, UTM 37S (primary); float-inclusive (sensitivity)

# %%
def bin5(s: pd.DataFrame, qmax: int) -> pd.DataFrame:
    s = s[s.Q <= qmax]
    x, y = TO_UTM.transform(s.lon.to_numpy(), s.lat.to_numpy())
    d = pd.DataFrame({"x": x, "y": y, "z": s.h.to_numpy()}, index=s.index)
    g = d.groupby(d.index.floor("5min"))
    m = g.median()
    return m[g.size() >= 150]


ds = {}
for name, s in sols.items():
    for tag, q in (("fixed", 1), ("float", 2)):
        b = bin5(s, q)
        for k in ("x", "y", "z"):
            ds[f"{k}_{name}_{tag}"] = b[k]
grid = pd.date_range("2021-12-31", "2022-01-07", freq="5min", inclusive="left")
kin = xr.Dataset({k: ("time", v.reindex(grid).to_numpy()) for k, v in ds.items()}, coords={"time": grid})
kin.attrs = {"title": "Kinematic GNSS, 5 min medians of 1 s epochs (UTM 37S, ellipsoidal height), UTC",
             "source": "ADS A20220506-004, Sugiyama, Minowa, Kondo & Aoki (2022), CC BY 4.0",
             "processing": "RTKLIB v2.5.1 rnx2rtkp kinematic, forward filter, base LGFX, GPS L1+L2, mask 15 deg, AR ratio 3"}
kin.to_netcdf(RESULTS / "kinematic_5min.nc")
print({k: int(np.isfinite(kin[k]).sum()) for k in kin.data_vars if k.endswith("fixed")})

# %% [markdown]
# ## Consistency check (stop rule, § 6.1)
#
# 15 min means of the kinematic positions vs the deposit's 15 min positions, each
# linearly detrended over the window: SD of the horizontal difference must be ≤ 5 cm.

# %%
rows = []
for name, mend in (("GNSS1", "gnss_GNSS1"), ("GNSS2", "gnss_GNSS2")):
    m = xr.open_dataset(ROOT / "data" / "clean" / f"{mend}.nc").to_dataframe()
    m.index = m.index.round("15min")
    m = m[~m.index.duplicated()]
    k = kin[[f"x_{name}_fixed", f"y_{name}_fixed", f"z_{name}_fixed"]].to_dataframe()
    k.columns = ["x", "y", "z"]
    k15 = k.groupby(k.index.floor("15min")).mean()
    j = k15.join(m[["x", "y", "z"]], rsuffix="_m", how="inner").dropna()
    t = (j.index - j.index[0]) / pd.Timedelta("1D")
    res = {}
    for c in ("x", "y", "z"):
        d = (j[c] - j[f"{c}_m"]).to_numpy()
        res[c] = d - np.polyval(np.polyfit(t, d, 1), t)
    sd_h = float(np.hypot(np.std(res["x"]), np.std(res["y"])))
    rows.append({"station": name, "n_15min": len(j), "sd_x_m": np.std(res["x"]), "sd_y_m": np.std(res["y"]),
                 "sd_z_m": np.std(res["z"]), "sd_horizontal_m": sd_h, "pass": sd_h <= 0.05})
check = pd.DataFrame(rows)
check.to_csv(RESULTS / "kinematic_consistency.csv", index=False)
print(check.round(4).to_string(index=False))
if not check["pass"].all():
    raise RuntimeError("Stop rule: kinematic positions disagree with the deposit (> 5 cm SD); review before use")
