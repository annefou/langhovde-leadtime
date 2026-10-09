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
# # 04 — Figures (exploratory, n = 1 event)
#
# `figures/main_result.png`: Period II series and the onset leads of basal pressure
# relative to glacier speed, under the three processings. The h = 6 h series are shown
# for display only. The lead panel uses every rule × bandwidth × reference window.

# %%
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

plt.style.use("seaborn-v0_8-whitegrid")
RESULTS, FIGURES = Path("../results"), Path("../figures")
FIGURES.mkdir(exist_ok=True)
s = xr.open_dataset(RESULTS / "series_15min.nc").sel(time=slice("2022-01-01", "2022-01-07"))
onsets = pd.read_csv(RESULTS / "onsets.csv", parse_dates=["onset"])
leads = pd.read_csv(RESULTS / "leads.csv")
delay = pd.read_csv(RESULTS / "robust_delay_check.csv")
H = 6
COL = {"centred": "#7f7f7f", "causal": "#d95f02", "causal-robust": "#1b9e77"}

# %%
fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(4, 2, width_ratios=[2.0, 1], hspace=0.15, wspace=0.55)
panels = [("tpos", "T⁺ (°C)"), ("level", "BH2201 level (m a.s.l.)"),
          ("uplift_GNSS1", "GNSS1 uplift (m)"), ("speed_GNSS1", "GNSS1 speed (m d⁻¹)")]
t = s.time.values
for i, (var, lab) in enumerate(panels):
    ax = fig.add_subplot(gs[i, 0])
    for proc in ["centred", "causal-robust"]:
        ax.plot(t, s[f"{var}_{proc}_h{H}"], color=COL[proc], lw=1.3, label=proc)
    if var == "tpos":
        ax2 = ax.twinx()
        ax2.fill_between(t, 0, s[f"rain_causal_h{H}"].clip(min=0), color="#7570b3", alpha=0.3, step="mid")
        ax2.set_ylabel("rain\n(mm h⁻¹)", color="#7570b3", fontsize=8)
        ax2.grid(False)
    o = onsets[(onsets.period == "II") & (onsets.series == var) & (onsets.h == H) & (onsets.ref == "R")
               & (onsets.rule == "CP")]
    for proc in ["centred", "causal-robust"]:
        oo = o[o.proc == proc].onset.dropna()
        for x in oo:
            ax.axvline(x, color=COL[proc], ls="--", lw=1)
    ax.axvspan(np.datetime64("2022-01-02"), np.datetime64("2022-01-06"), color="0.9", zorder=-1)
    ax.set_ylabel(lab, fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d Jan"))
    if i < 3:
        ax.set_xticklabels([])
    if i == 0:
        ax.legend(loc="upper right", fontsize=8, frameon=True)
        ax.set_title(f"Period II (shaded), h = {H} h; dashed = change-point onset", fontsize=10)

ax = fig.add_subplot(gs[:, 1])
pairs = [("level", "speed_GNSS1"), ("level", "accel_GNSS1"), ("level", "speed_GNSS2"),
         ("uplift_GNSS1", "speed_GNSS1")]
rules = ["T2", "T3", "T5", "CP"]
y = 0
yt, yl = [], []
for b, e in pairs:
    for proc in ["centred", "causal", "causal-robust"]:
        v = leads[(leads.basal_or_forcing == b) & (leads.endpoint == e) & (leads.proc == proc)
                  & leads.rule.isin(rules)].lead_h.dropna()
        ax.scatter(v, np.full(v.size, y) + np.random.default_rng(0).uniform(-0.25, 0.25, v.size),
                   s=10, color=COL[proc], alpha=0.7)
        if v.size:
            ax.plot([v.median()] * 2, [y - 0.35, y + 0.35], color="k", lw=1.5)
        yt.append(y)
        yl.append(f"{b.replace('_GNSS1', '')} → {e}\n{proc} (n={v.size}/32)")
        y += 1
    y += 0.6
lo, hi = delay.loc[delay.h >= 3, "cp_delay_median_h_causal-robust"].agg(["min", "max"])
ax.axvspan(lo, hi, color=COL["causal-robust"], alpha=0.12)
ax.axvspan(-1, 1, color="0.85")
ax.axvline(0, color="k", lw=0.8)
ax.set_yticks(yt, yl, fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel("lead (h): > 0 = basal series changes first", fontsize=9)
ax.set_title("Onset leads, Period II (n = 1 event)\ngrey ±1 h; green = causal-robust detection delay",
             fontsize=9)
fig.suptitle("Langhovde, 2–6 Jan 2022: does the bed signal come first? One event; exploratory",
             fontsize=11)
fig.savefig(FIGURES / "main_result.png", dpi=150, bbox_inches="tight")
fig.savefig(FIGURES / "main_result.pdf", bbox_inches="tight")
plt.show()
