"""Guards for notebooks/leadtime.py (ANALYSIS_PLAN.md § 1)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))
from leadtime import (backward_rate, local_linear, onset_changepoint,  # noqa: E402
                      onset_threshold, smooth_on_grid, xcorr)

PRIOR_GNSS = ROOT / "data" / "raw" / "prior_v1.0.2" / "notebooks"
rng = np.random.default_rng(0)


def test_causal_ignores_future_samples():
    t = np.arange(0, 10, 0.01)
    x = np.sin(t) + rng.normal(0, 0.1, t.size)
    ti = np.array([3.0, 5.0, 7.0])
    base = local_linear(t, x, ti, 0.5, causal=True)
    x2 = x.copy()
    x2[t > 5.0] += 100.0
    changed = local_linear(t, x2, ti, 0.5, causal=True)
    assert changed[0] == pytest.approx(base[0], abs=0)
    assert changed[1] == pytest.approx(base[1], abs=0)
    assert abs(changed[2] - base[2]) > 1


def test_causal_local_linear_tracks_a_ramp_without_lag():
    t = np.arange(0, 10, 0.01)
    x = 2.0 * t + 1.0
    ti = np.array([4.0, 6.0])
    assert local_linear(t, x, ti, 1.0, causal=True) == pytest.approx(2 * ti + 1, abs=1e-9)
    # The local-constant one-sided fit lags a ramp by h·√(2/π) (amendment A1).
    lagged = local_linear(t, x, ti, 1.0, causal=True, degree=0)
    assert (2 * ti + 1 - lagged) / 2 == pytest.approx(np.sqrt(2 / np.pi), abs=0.02)


def test_min_history():
    t = np.arange(0, 1, 0.01)
    out = local_linear(t, t, np.array([0.05, 0.5]), 0.1, causal=True, min_history=0.1)
    assert np.isnan(out[0]) and np.isfinite(out[1])


@pytest.mark.skipif(not (PRIOR_GNSS / "gnss.py").exists(), reason="run 01_data_download first")
def test_centred_branch_matches_prior_local_regression():
    sys.path.insert(0, str(PRIOR_GNSS))
    from gnss import local_regression
    t = np.sort(rng.uniform(0, 5, 400))
    x = np.cos(3 * t) + rng.normal(0, 0.2, t.size)
    ti = np.linspace(0.5, 4.5, 50)
    assert local_linear(t, x, ti, 0.2, causal=False) == pytest.approx(local_regression(t, x, ti, 0.2), abs=1e-9)


def test_smooth_on_grid_respects_segments():
    time = pd.date_range("2022-01-01", periods=200, freq="15min")
    time = time.delete(range(80, 140))  # 15 h gap
    grid = pd.date_range("2022-01-01", periods=200, freq="15min")
    out = smooth_on_grid(time.values, np.ones(time.size), grid, 1.0, causal=True)
    assert np.isnan(out[85]) and np.isnan(out[140]) and np.isfinite(out[150])


def test_backward_rate():
    x = np.arange(10.0)
    r = backward_rate(x, step_h=0.25)
    assert np.isnan(r[0]) and r[1:] == pytest.approx(96.0)


def _event_series(onset: str) -> pd.Series:
    idx = pd.date_range("2022-01-01", "2022-01-06", freq="15min")
    t = (idx - pd.Timestamp(onset)) / pd.Timedelta("1h")
    return pd.Series(np.clip(t, 0, None) * 0.5 + rng.normal(0, 0.05, idx.size), index=idx)


def test_onsets_find_a_known_hinge():
    s = _event_series("2022-01-03 06:00")
    ref, search = ("2022-01-01", "2022-01-02"), ("2022-01-02", "2022-01-06")
    thr = onset_threshold(s, ref, search, k=3)
    assert pd.Timestamp("2022-01-03 06:00") <= thr <= pd.Timestamp("2022-01-03 07:00")
    tau, lo, hi = onset_changepoint(s, "2022-01-01", search)
    assert abs(tau - pd.Timestamp("2022-01-03 06:00")) <= pd.Timedelta("30min")
    assert lo <= tau <= hi


def test_xcorr_sign_convention():
    idx = pd.date_range("2022-01-01", periods=2000, freq="15min")
    x = pd.Series(rng.normal(size=idx.size), index=idx).rolling(8, min_periods=1).mean()
    y = x.shift(12)  # y is x delayed by 3 h: x leads y
    r = xcorr(x, y, (str(idx[0]), str(idx[-1])), 24)
    assert r.idxmax() == 12
