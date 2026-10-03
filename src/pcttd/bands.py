"""Band-structure post-processing: gaps, guided-mode counting, single-mode windows, slow-light metrics.

Pure NumPy (no MPB) so every definition used in the Phase 0.5 geometry-selection criteria is unit-testable.

Units: w = omega*a/(2 pi c), k = k_phys*a/(2 pi); the group velocity in units of c is v = dw/dk,
and the group index is n_g = 1/v.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# --------------------------------------------------------------------------- bulk gaps
@dataclass(frozen=True)
class Gap:
    n: int            # gap lies between band n and band n+1 (1-based)
    w_lo: float
    w_hi: float

    @property
    def w_mid(self) -> float:
        return 0.5 * (self.w_lo + self.w_hi)

    @property
    def width(self) -> float:
        return self.w_hi - self.w_lo

    @property
    def rel(self) -> float:
        return self.width / self.w_mid


def find_gaps(freqs: np.ndarray, n_max: int | None = None) -> list[Gap]:
    """Complete gaps from freqs[k, band] (all k of the BZ sample). Only strictly positive gaps are returned."""
    freqs = np.asarray(freqs)
    nb = freqs.shape[1]
    top = nb - 1 if n_max is None else min(n_max, nb - 1)
    gaps = []
    for n in range(1, top + 1):
        lo = freqs[:, n - 1].max()
        hi = freqs[:, n].min()
        if hi > lo:
            gaps.append(Gap(n, float(lo), float(hi)))
    return gaps


# --------------------------------------------------------------------------- Hermite resampling
def hermite_resample(k: np.ndarray, w: np.ndarray, dw: np.ndarray, per_interval: int = 24):
    """Cubic-Hermite resample of w(k) given values and slopes on a coarse grid.

    Returns (kf, wf, dwf) on a fine grid (analytic derivative of the Hermite cubic).
    """
    k, w, dw = (np.asarray(a, float) for a in (k, w, dw))
    kf, wf, df = [], [], []
    t = np.linspace(0.0, 1.0, per_interval, endpoint=False)
    h00 = 2 * t**3 - 3 * t**2 + 1
    h10 = t**3 - 2 * t**2 + t
    h01 = -2 * t**3 + 3 * t**2
    h11 = t**3 - t**2
    d00 = 6 * t**2 - 6 * t
    d10 = 3 * t**2 - 4 * t + 1
    d01 = -6 * t**2 + 6 * t
    d11 = 3 * t**2 - 2 * t
    for i in range(len(k) - 1):
        h = k[i + 1] - k[i]
        kf.append(k[i] + t * h)
        wf.append(h00 * w[i] + h10 * h * dw[i] + h01 * w[i + 1] + h11 * h * dw[i + 1])
        df.append((d00 * w[i] + d01 * w[i + 1]) / h + d10 * dw[i] + d11 * dw[i + 1])
    kf.append([k[-1]])
    wf.append([w[-1]])
    df.append([dw[-1]])
    return np.concatenate(kf), np.concatenate(wf), np.concatenate(df)


# --------------------------------------------------------------------------- crossings
@dataclass
class Crossings:
    """All band crossings of the level grid. One record per (level, band-branch crossing)."""
    levels: np.ndarray
    level_idx: np.ndarray
    vg: np.ndarray
    parity: np.ndarray    # +1 even-Phi, -1 odd-Phi
    strip: np.ndarray
    edge: np.ndarray
    k: np.ndarray

    def counts(self, mask: np.ndarray | None = None) -> np.ndarray:
        idx = self.level_idx if mask is None else self.level_idx[mask]
        return np.bincount(idx, minlength=len(self.levels))


def class_branches(freq, vg, parity, strip, edge):
    """Split bands at each k into parity classes and sort within each class by frequency.

    Sorting within a parity class removes the spurious kinks that sorting across classes would
    create at (symmetry-allowed) even/odd crossings. Returns a list of dicts, one per (class, branch).
    Only bands below the lowest 'highest computed band' are kept so the branch count is k-independent.
    """
    freq, vg, parity, strip, edge = (np.asarray(a) for a in (freq, vg, parity, strip, edge))
    nk, nb = freq.shape
    wcut = freq[:, -1].min()
    out = []
    for p in (+1, -1):
        orders, counts = [], []
        for ik in range(nk):
            sel = np.where((parity[ik] == p) & (freq[ik] < wcut))[0]
            sel = sel[np.argsort(freq[ik, sel])]
            orders.append(sel)
            counts.append(len(sel))
        m = min(counts)
        for b in range(m):
            rows = np.arange(nk)
            cols = np.array([orders[ik][b] for ik in range(nk)])
            out.append(dict(parity=p, w=freq[rows, cols], vg=vg[rows, cols],
                            strip=strip[rows, cols], edge=edge[rows, cols]))
    return out


def find_crossings(k, branches, levels, per_interval: int = 24) -> Crossings:
    levels = np.asarray(levels, float)
    li, vgs, par, st, ed, kk = [], [], [], [], [], []
    for br in branches:
        kf, wf, df = hermite_resample(k, br["w"], br["vg"], per_interval)
        sf = np.interp(kf, k, br["strip"])
        ef = np.interp(kf, k, br["edge"])
        lo = np.minimum(wf[:-1], wf[1:])
        hi = np.maximum(wf[:-1], wf[1:])
        i0 = np.searchsorted(levels, lo, side="left")
        i1 = np.searchsorted(levels, hi, side="left")  # level < hi, so a shared fine-grid endpoint counts once
        for j in np.nonzero(i1 > i0)[0]:
            for L in range(i0[j], i1[j]):
                dw = wf[j + 1] - wf[j]
                if dw == 0:
                    continue
                t = (levels[L] - wf[j]) / dw
                li.append(L)
                vgs.append(df[j] + t * (df[j + 1] - df[j]))
                par.append(br["parity"])
                st.append(sf[j] + t * (sf[j + 1] - sf[j]))
                ed.append(ef[j] + t * (ef[j + 1] - ef[j]))
                kk.append(kf[j] + t * (kf[j + 1] - kf[j]))
    return Crossings(levels, np.array(li, int), np.array(vgs), np.array(par, int),
                     np.array(st), np.array(ed), np.array(kk))


# --------------------------------------------------------------------------- windows
def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous True runs as (start, stop) index pairs, stop exclusive."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(int))
    return list(zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]))


def single_mode_masks(c: Crossings):
    """Per-level masks (Revision 1 definitions; see docs/phase0p5_criteria.md).

    A band crossing of the level over k in [0, 0.5] is exactly one forward-propagating guided mode,
    whatever the sign of dw/dk in the folded zone: its time-reversed partner (-k, -v_g) carries power the
    other way, and a negative slope on [0, 0.5] is a positive-power-flow mode with true Bloch wavevector in
    (0.5, 1). Hence the number of forward modes at level w equals the number of crossings.

    strict_any  : exactly one forward mode (either parity).
    strict_even : exactly one forward mode in total, and it is even-Phi.
    even_only   : exactly one even-Phi forward mode (odd modes may coexist).
    """
    tot = c.counts()
    ev = c.counts(c.parity == +1)
    strict_any = tot == 1
    strict_even = strict_any & (ev == 1)
    even_only = ev == 1
    return strict_any, strict_even, even_only


def per_level_single(c: Crossings, rec_mask: np.ndarray):
    """Per-level v_g, strip and edge energy fractions taken from the crossing records selected by ``rec_mask``
    (boolean per record; at most one selected record per level is expected). nan where no record."""
    n = len(c.levels)
    vg = np.full(n, np.nan)
    strip = np.full(n, np.nan)
    edge = np.full(n, np.nan)
    sel = np.nonzero(rec_mask)[0]
    vg[c.level_idx[sel]] = c.vg[sel]
    strip[c.level_idx[sel]] = c.strip[sel]
    edge[c.level_idx[sel]] = c.edge[sel]
    return vg, strip, edge


def frac_width(levels: np.ndarray, start: int, stop: int) -> float:
    lo, hi = levels[start], levels[stop - 1]
    return float((hi - lo) / (0.5 * (hi + lo)))


# --------------------------------------------------------------------------- slow light
def best_flat_window(w: np.ndarray, ng: np.ndarray, n_min: float, tol: float = 0.10,
                     max_pts: int = 300):
    """Largest group-index-bandwidth product over sub-intervals with n_g within +-tol of their mean.

    Returns dict(gbp, n_mean, i, j, w_lo, w_hi) or None if no interval qualifies. ``w`` must be uniformly
    spaced and increasing (the level grid). Mean is the frequency average.
    GBP = n_mean * (w_hi - w_lo) / w_center.
    """
    w = np.asarray(w, float)
    ng = np.asarray(ng, float)
    if len(w) > max_pts:
        idx = np.unique(np.linspace(0, len(w) - 1, max_pts).round().astype(int))
        w, ng = w[idx], ng[idx]
    n = len(w)
    best = None
    cs = np.concatenate([[0.0], np.cumsum(ng)])
    for i in range(n - 1):
        j = np.arange(i + 1, n)
        cnt = j - i + 1
        mean = (cs[j + 1] - cs[i]) / cnt
        mx = np.maximum.accumulate(ng[i + 1:])
        mx = np.maximum(mx, ng[i])
        mn = np.minimum.accumulate(ng[i + 1:])
        mn = np.minimum(mn, ng[i])
        ok = (mx <= (1 + tol) * mean) & (mn >= (1 - tol) * mean) & (mean >= n_min)
        if not ok.any():
            continue
        wc = 0.5 * (w[j] + w[i])
        gbp = mean * (w[j] - w[i]) / wc
        gbp = np.where(ok, gbp, -np.inf)
        jj = int(np.argmax(gbp))
        if np.isfinite(gbp[jj]) and (best is None or gbp[jj] > best["gbp"]):
            best = dict(gbp=float(gbp[jj]), n_mean=float(mean[jj]), w_lo=float(w[i]),
                        w_hi=float(w[j[jj]]), frac_bw=float((w[j[jj]] - w[i]) / wc[jj]))
    return best
