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
# # 02 — Data clean (Step 0)
#
# Runs the prior chain's own parsing code rather than rewriting it: the cells of
# v1.0.2 `notebooks/02_data_clean.py` (Zenodo
# [doi:10.5281/zenodo.23257925](https://doi.org/10.5281/zenodo.23257925)) that parse
# the Fig. 3 inputs of the Mendeley deposit. The later cells of that notebook parse
# JMA and NOAA records, which this analysis does not use, and are not run.
#
# Output in `data/clean/`: `gnss_GNSS1.nc`, `gnss_GNSS2.nc`, `pressure_BH22*.nc`,
# `tide.nc`, `aws.nc` (NetCDF, CF-1.8; all times UTC as established in v1.0.2).

# %%
from pathlib import Path

PRIOR = Path("../data/raw/prior_v1.0.2")
src = (PRIOR / "notebooks" / "02_data_clean.py").read_text()
STOP = "# %% [markdown]\n# ## Syowa daily temperature"
if STOP not in src:
    raise RuntimeError("v1.0.2 02_data_clean.py layout changed; cannot locate the Fig. 3 cells")
fig3_cells = src.split(STOP)[0]

# The prior notebook resolves paths relative to its own notebooks/ directory,
# which is the same layout as here (../data/raw/mendeley, ../data/clean).
exec(compile(fig3_cells, str(PRIOR / "notebooks" / "02_data_clean.py"), "exec"), {"__name__": "__prior02__"})

# %%
expected = ["gnss_GNSS1.nc", "gnss_GNSS2.nc", "pressure_BH2201.nc", "pressure_BH2202.nc",
            "pressure_BH2203.nc", "tide.nc", "aws.nc"]
missing = [f for f in expected if not (Path("../data/clean") / f).exists()]
if missing:
    raise RuntimeError(f"not produced: {missing}")
print("clean files:", expected)
