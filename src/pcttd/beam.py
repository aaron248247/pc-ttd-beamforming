"""Phase 5 beam metrics built on ``pcttd.array`` (pure NumPy): lobe finding on the full visible region, local HPBW of the design lobe,
phase-error-to-pointing link and per-case analysis over a frequency set.

The design lobe is the principal maximum inside |sin(theta) - sin(theta_0)| <= window (``array.pointing_angle``); every other local maximum of the
power over -90..+90 degrees is a sidelobe or a grating lobe.  All power values are |AF|^2 / N^2 (1 for an ideal equal-amplitude array at its peak).
"""
from __future__ import annotations

import numpy as np

from . import array as A

S_GRID = 20001                       # full-region grid in sin(theta): ds = 1e-4


def full_power(H_f, f, d, n=S_GRID, c=A.C):
    """(sin_theta grid, power |AF|^2/N^2) over the whole visible region for channel weights H_f[n] at frequency f."""
    s = np.linspace(-1.0, 1.0, n)
    AF = A.array_factor(np.asarray(H_f, complex)[:, None], [f], d, s, c)[0]
    return s, np.abs(AF) ** 2 / len(H_f) ** 2


def local_maxima(p):
    """Indices of local maxima of p (plateaus count once at their first sample); the two end samples count if higher than their neighbor."""
    p = np.asarray(p, float)
    idx = [i for i in range(1, len(p) - 1) if p[i] > p[i - 1] and p[i] >= p[i + 1]]
    if p[0] > p[1]:
        idx.insert(0, 0)
    if p[-1] > p[-2]:
        idx.append(len(p) - 1)
    return np.array(idx, int)


def lobe_report(H_f, f, d, theta_peak_deg, n=S_GRID):
    """Lobe analysis at one frequency.  The design lobe is the local maximum nearest sin(theta_peak_deg).
    Returns dict: gain_design, hpbw_deg, strongest_db (strongest other local maximum relative to the design peak, dB; -inf if none),
    strongest_theta_deg, n_within_3db (interior local maxima within 3 dB of the design peak, design lobe included; the +-90 deg endpoints are not counted), edge_db (power at -90 and +90 deg relative to the design peak), lobes (theta_deg, level_db arrays)."""
    s, p = full_power(H_f, f, d, n)
    mx = local_maxima(p)
    sp = np.sin(np.deg2rad(theta_peak_deg))
    k = mx[np.argmin(np.abs(s[mx] - sp))]
    pk = p[k]
    other = mx[mx != k]
    lev = 10 * np.log10(np.maximum(p[other], 1e-300) / pk)
    j = int(np.argmax(lev)) if len(other) else None
    return dict(gain_design=float(pk), hpbw_deg=_hpbw(s, p, k),
                strongest_db=float(lev[j]) if j is not None else -np.inf,
                strongest_theta_deg=float(np.rad2deg(np.arcsin(s[other[j]]))) if j is not None else np.nan,
                n_within_3db=int(np.sum(p[mx[(mx > 0) & (mx < len(p) - 1)]] >= pk / 2)),
                edge_db=(float(10 * np.log10(p[0] / pk)), float(10 * np.log10(p[-1] / pk))),
                lobes=(np.rad2deg(np.arcsin(s[mx])), 10 * np.log10(np.maximum(p[mx], 1e-300) / pk)))


def _crossing(s, p, i_from, i_to, level):
    step = 1 if i_to > i_from else -1
    j = i_from
    while j != i_to and p[j] > level:
        j += step
    if j == i_to and p[j] > level:
        return np.nan
    j0 = j - step                                       # last sample above the level
    return s[j0] + (s[j] - s[j0]) * (p[j0] - level) / (p[j0] - p[j])


def _hpbw(s, p, k):
    lo = _crossing(s, p, k, 0, p[k] / 2)
    hi = _crossing(s, p, k, len(p) - 1, p[k] / 2)
    if np.isnan(lo) or np.isnan(hi):
        return float("nan")
    return float(np.rad2deg(np.arcsin(hi) - np.arcsin(lo)))


def analyze_case(H, freqs, d, theta0_deg, window, n_grid=40001):
    """Per-frequency design-lobe metrics for H[n, f] on ``freqs``.  Returns dict of arrays over freqs."""
    freqs = np.asarray(freqs, float)
    th, gain, _ = A.pointing_angle(H, freqs, d, theta0_deg, window=window, n_grid=n_grid)
    out = dict(theta_peak=th, gain=gain, hpbw=np.empty(len(freqs)), strongest_db=np.empty(len(freqs)),
               strongest_theta=np.empty(len(freqs)), n_3db=np.empty(len(freqs), int), edge_db=np.empty((len(freqs), 2)))
    for i, f in enumerate(freqs):
        r = lobe_report(H[:, i], f, d, th[i])
        out["hpbw"][i], out["strongest_db"][i], out["strongest_theta"][i], out["n_3db"][i] = r["hpbw_deg"], r["strongest_db"], r["strongest_theta_deg"], r["n_within_3db"]
        out["edge_db"][i] = r["edge_db"]
    return out


def grating_lobe_angles(theta_deg, freqs, d, c=A.C):
    """Analytic positions (deg) of the visible replicas of a lobe at theta_deg: sin = sin(theta) + q c/(f d), q != 0, |sin| <= 1.  List per frequency."""
    s0 = np.sin(np.deg2rad(np.broadcast_to(theta_deg, np.shape(freqs))))
    res = []
    for s, f in zip(s0, freqs):
        q = np.arange(-8, 9)
        q = q[q != 0]
        ss = s + q * c / (f * d)
        res.append(np.sort(np.rad2deg(np.arcsin(ss[np.abs(ss) <= 1]))))
    return res


def phase_slope_pointing(eps, freqs, d, theta0_deg, c=A.C):
    """Linear-model pointing shift (deg) caused by residual phase errors eps[n, f] (rad) relative to the ideal weights:
    least-squares slope g(f) of eps vs element index; Delta(sin theta) = g c / (2 pi f d); Delta theta = Delta(sin theta) / cos(theta_0)."""
    eps = np.asarray(eps, float)
    m = np.arange(eps.shape[0], dtype=float)
    mc = m - m.mean()
    g = (mc[:, None] * eps).sum(0) / (mc ** 2).sum()
    return np.rad2deg(g * c / (2 * np.pi * np.asarray(freqs) * d) / np.cos(np.deg2rad(theta0_deg))), g
