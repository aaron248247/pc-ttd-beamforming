"""True-time-delay figures of merit (Phase 3+), written from scratch after the crash-recovery stubs were deleted.

Pure NumPy; no solver dependency.  Conventions (a = c = 1, f = a/lambda, omega = 2 pi f):
    tau(f)        per-period group delay [a/c] = d theta_pp / d omega  (equals the group index n_g in these units)
    GDD(f)        per-period group-delay dispersion d tau / d omega   [(a/c)^2]
    ripple        (max - min) / median of tau over an interval          (dimensionless)
    frac. BW      (f_high - f_low) / f_center, f_center = (f_low + f_high) / 2
    GBP           fractional bandwidth x median n_g (the declared group-delay-bandwidth product)
Phase extraction, unwrapping and the cascade de-embedding live in ``guided`` / ``w1analysis`` and are not duplicated here.
"""
from __future__ import annotations

import numpy as np

from . import extract as X


def omega(f):
    return X.omega(np.asarray(f, float))


def fractional_bandwidth(f_low, f_high):
    """(f_high - f_low) / f_center with f_center the arithmetic midpoint."""
    f_low, f_high = float(f_low), float(f_high)
    if f_high < f_low:
        raise ValueError("f_high < f_low")
    return (f_high - f_low) / (0.5 * (f_low + f_high))


def gdd(f, tau):
    """d tau / d omega by second-order central differences (one-sided at the ends), omega = 2 pi f."""
    return np.gradient(np.asarray(tau, float), omega(f))


def median_delay(tau, mask=None):
    t = np.asarray(tau, float)
    if mask is not None:
        t = t[np.asarray(mask, bool)]
    return float(np.nanmedian(t))


def delay_ripple(tau, mask=None):
    """(max - min) / median over the (masked) samples."""
    t = np.asarray(tau, float)
    if mask is not None:
        t = t[np.asarray(mask, bool)]
    med = np.nanmedian(t)
    return float((np.nanmax(t) - np.nanmin(t)) / med)


def max_relative_deviation(tau, mask=None):
    """max |tau / median - 1| over the (masked) samples (the quantity bounded by the +-tol flatness rule)."""
    t = np.asarray(tau, float)
    if mask is not None:
        t = t[np.asarray(mask, bool)]
    med = np.nanmedian(t)
    return float(np.nanmax(np.abs(t / med - 1.0)))


def contiguous_runs(mask):
    """Index pairs (i0, i1) inclusive of the maximal runs of True in ``mask``."""
    m = np.asarray(mask, bool)
    runs, start = [], None
    for i, v in enumerate(m):
        if v and start is None:
            start = i
        if not v and start is not None:
            runs.append((start, i - 1)); start = None
    if start is not None:
        runs.append((start, len(m) - 1))
    return runs


def usable_band(f, tau, ok, tol=0.10, min_samples=2):
    """Widest contiguous frequency interval on which every sample is ``ok`` (e.g. T >= 0.5 at every length, inside the allowed domain,
    finite tau) AND tau stays within +-tol of its own median over the interval.

    The median depends on the interval, so all windows are tested (brute force; the flatness condition is not monotone in the window).
    Ties in width are broken by the smaller ripple.  Returns None if no window with ``min_samples`` samples qualifies, else a dict with
    indices, f_low, f_high, f_center, frac_bw, median_tau, ripple, max_dev."""
    f = np.asarray(f, float); tau = np.asarray(tau, float)
    ok = np.asarray(ok, bool) & np.isfinite(tau)
    best = None
    for i0, i1 in contiguous_runs(ok):
        for a in range(i0, i1 + 1):
            for b in range(a + min_samples - 1, i1 + 1):
                seg = tau[a:b + 1]
                med = np.median(seg)
                if med <= 0:
                    continue
                dev = float(np.max(np.abs(seg / med - 1.0)))
                if dev > tol:
                    continue
                width = f[b] - f[a]
                rip = float((seg.max() - seg.min()) / med)
                key = (width, -rip)
                if best is None or key > best[0]:
                    best = (key, a, b, med, rip, dev)
    if best is None:
        return None
    _, a, b, med, rip, dev = best
    return dict(i_low=int(a), i_high=int(b), f_low=float(f[a]), f_high=float(f[b]), f_center=float(0.5 * (f[a] + f[b])),
                frac_bw=fractional_bandwidth(f[a], f[b]), median_tau=float(med), ripple=rip, max_dev=dev, tol=float(tol))


def delay_bandwidth_product(frac_bw, median_ng):
    """Declared group-delay-bandwidth product: fractional bandwidth x median group index (= median per-period delay in a/c)."""
    return float(frac_bw) * float(median_ng)


def required_length_periods(tau_target, tau_per_period):
    """Number of periods (continuous) needed to reach a total delay ``tau_target`` [a/c] with per-period delay ``tau_per_period`` [a/c]."""
    if tau_per_period <= 0:
        raise ValueError("per-period delay must be positive")
    return float(tau_target) / float(tau_per_period)


def target_delay_from_cycles(n_cycles, f_center):
    """Delay of ``n_cycles`` carrier periods at ``f_center`` [a/c]: n_cycles / f_center."""
    return float(n_cycles) / float(f_center)


def common_band(bands):
    """Intersection of (f_low, f_high) intervals; None entries or an empty intersection give None."""
    if any(b is None for b in bands):
        return None
    lo = max(b[0] for b in bands); hi = min(b[1] for b in bands)
    return (float(lo), float(hi)) if hi >= lo else None


def interp_on_grid(f_new, f, y):
    """Linear interpolation of a real sample set onto ``f_new`` (NaN outside the sampled range; no extrapolation)."""
    return np.interp(np.asarray(f_new, float), np.asarray(f, float), np.asarray(y, float), left=np.nan, right=np.nan)
