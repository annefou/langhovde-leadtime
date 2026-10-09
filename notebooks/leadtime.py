"""Causal (one-sided) smoothing, onset rules and lead-lag tools for ANALYSIS_PLAN.md.

The centred smoother is the prior chain's `local_regression` (v1.0.2 `gnss.py`,
imported unmodified). The causal smoother here is the same Gaussian-kernel local
linear regression with zero weight on samples later than the evaluation time.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

_DAY = np.timedelta64(1, "D")


def _days(t: np.ndarray, t0: np.datetime64) -> np.ndarray:
    return (np.asarray(t, dtype="datetime64[ns]") - t0) / _DAY


def local_linear(t: np.ndarray, x: np.ndarray, ti: np.ndarray, h: float, causal: bool,
                 degree: int = 1, min_history: float = 0.0) -> np.ndarray:
    """Gaussian-kernel (σ = h) local polynomial fit of x(t) at ti; degree 1 or 0.

    `causal=True` gives zero weight to samples later than ti and needs samples
    covering at least `min_history` before ti. Units of t, ti, h, min_history agree.
    Samples beyond 8σ are dropped (weight < 2e-14), which only speeds things up.
    """
    t, x, ti = (np.asarray(a, dtype=float) for a in (t, x, ti))
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    order = np.argsort(t)
    t, x = t[order], x[order]
    lo = np.searchsorted(t, ti - 8 * h, side="left")
    hi = np.searchsorted(t, ti if causal else ti + 8 * h, side="right")
    out = np.full(ti.size, np.nan)
    for k, t0 in enumerate(ti):
        tt, xx = t[lo[k]:hi[k]], x[lo[k]:hi[k]]
        if tt.size < 2 or (causal and t0 - tt[0] < min_history):
            continue
        d = tt - t0
        w = np.exp(-0.5 * (d / h) ** 2)
        if degree == 0:
            out[k] = (w * xx).sum() / w.sum()
            continue
        sw, swd, swdd = w.sum(), (w * d).sum(), (w * d * d).sum()
        det = sw * swdd - swd ** 2
        if det > 1e-300:
            out[k] = (swdd * (w * xx).sum() - swd * (w * d * xx).sum()) / det
    return out


def segments(time: np.ndarray, max_gap_h: float) -> list[tuple[np.datetime64, np.datetime64]]:
    """(start, end) of runs with no sampling gap longer than max_gap_h."""
    time = np.asarray(time, dtype="datetime64[ns]")
    dt_h = np.diff(time) / np.timedelta64(1, "h")
    breaks = np.flatnonzero(dt_h > max_gap_h) + 1
    edges = np.concatenate([[0], breaks, [time.size]])
    return [(time[a], time[b - 1]) for a, b in zip(edges[:-1], edges[1:])]


def smooth_on_grid(time: np.ndarray, values: np.ndarray, grid: pd.DatetimeIndex, h_hours: float,
                   causal: bool, max_gap_h: float = 12.0, min_history_h: float = 1.0,
                   degree: int = 1) -> np.ndarray:
    """Smoothed series on `grid`, segment by segment; empty outside segments.

    Each segment is fitted from its own samples only. Causal values need
    `min_history_h` of the current segment behind them.
    """
    time = np.asarray(time, dtype="datetime64[ns]")
    values = np.asarray(values, dtype=float)
    g = grid.values.astype("datetime64[ns]")
    t0 = g[0]
    out = np.full(g.size, np.nan)
    for a, b in segments(time, max_gap_h):
        sel = (time >= a) & (time <= b)
        gi = (g >= a) & (g <= b)
        if not gi.any():
            continue
        out[gi] = local_linear(_days(time[sel], t0), values[sel], _days(g[gi], t0), h_hours / 24,
                               causal, degree, min_history_h / 24)
    return out


def backward_rate(x: np.ndarray, step_h: float, lag_steps: int = 1) -> np.ndarray:
    """(x(t) − x(t − lag)) / lag per day, using only the past."""
    out = np.full(x.size, np.nan)
    out[lag_steps:] = (x[lag_steps:] - x[:-lag_steps]) / (lag_steps * step_h / 24)
    return out


# ---------- onsets ----------

def onset_threshold(s: pd.Series, ref: tuple[str, str], search: tuple[str, str], k: float,
                    persist: int = 5) -> pd.Timestamp | None:
    """First time in `search` where s > μ + kσ (from `ref`) for `persist` consecutive points."""
    r = s[ref[0]:ref[1]].dropna()
    if r.size < 10:
        return None
    level = r.mean() + k * r.std(ddof=1)
    w = s[search[0]:search[1]]
    above = (w > level).to_numpy()
    run = 0
    for i, flag in enumerate(above):
        run = run + 1 if flag else 0
        if run == persist:
            return w.index[i - persist + 1]
    return None


def onset_changepoint(s: pd.Series, start: str, search: tuple[str, str], min_side: int = 4
                      ) -> tuple[pd.Timestamp, pd.Timestamp, pd.Timestamp] | None:
    """Hinge fit (flat, then linear) on [start, argmax in search]: (τ̂, low, high).

    The interval is all τ with SSE ≤ SSE_min·(1 + F₀.₉₅(1, n−3)/(n−3)); it ignores
    autocorrelation and is too narrow.
    """
    w = s[search[0]:search[1]].dropna()
    if w.empty:
        return None
    y = s[start:w.idxmax()].dropna()
    n = y.size
    if n < 2 * min_side + 3:
        return None
    t = (y.index - y.index[0]) / pd.Timedelta("1h")
    t = t.to_numpy(float)
    yv = y.to_numpy(float)
    taus, sse = [], []
    for i in range(min_side, n - min_side):
        X = np.column_stack([np.ones(n), np.clip(t - t[i], 0, None)])
        beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
        if beta[1] <= 0:
            continue
        taus.append(i)
        sse.append(((yv - X @ beta) ** 2).sum())
    if not taus:
        return None
    sse = np.array(sse)
    best = int(np.argmin(sse))
    crit = sse[best] * (1 + stats.f.ppf(0.95, 1, n - 3) / (n - 3))
    inside = [taus[j] for j in np.flatnonzero(sse <= crit)]
    return y.index[taus[best]], y.index[min(inside)], y.index[max(inside)]


# ---------- lead-lag ----------

def xcorr(x: pd.Series, y: pd.Series, window: tuple[str, str], max_lag_steps: int
          ) -> pd.Series:
    """Pearson r(ℓ) = corr(x(t), y(t + ℓ)) in `window`, each detrended; ℓ in steps.

    Positive ℓ means x leads y. Pairs with a missing value are dropped per lag.
    """
    def detrend(s: pd.Series) -> pd.Series:
        s = s[window[0]:window[1]]
        ok = s.notna()
        if ok.sum() < 3:
            return s
        tt = np.arange(s.size, dtype=float)
        p = np.polyfit(tt[ok], s[ok], 1)
        return s - np.polyval(p, tt)

    xd, yd = detrend(x), detrend(y).reindex(x[window[0]:window[1]].index)
    xv, yv = xd.to_numpy(float), yd.to_numpy(float)
    out = {}
    for lag in range(-max_lag_steps, max_lag_steps + 1):
        if lag >= 0:
            a, b = xv[: xv.size - lag], yv[lag:]
        else:
            a, b = xv[-lag:], yv[: yv.size + lag]
        ok = np.isfinite(a) & np.isfinite(b)
        out[lag] = np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 10 else np.nan
    return pd.Series(out)
