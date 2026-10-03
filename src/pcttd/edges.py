"""Stop-band edges of a periodic structure from the attenuation per period of finite sections.

Why not read the edge off T(f) of one finite section? The transmission of a finite stack near a band edge is dominated by
interface (impedance) mismatch and Fabry-Perot fringes, and a fixed threshold such as -3 dB drifts with length. The attenuation per
period, kappa(f) = -d ln T / dN, is a property of the periodic medium alone: interfaces contribute a length-independent factor that
the slope removes. Inside a gap kappa = 2 Im(K a) > 0; in a pass band it vanishes (T oscillates with N but does not decay).

Near a parabolic band edge, Im(K) ~ sqrt(|f - f_edge|), so kappa^2 is linear in f. The edge is the zero of a straight-line fit of kappa^2 vs f
over a window of moderately small kappa. The method's bias (curvature, window choice) is measured by applying the identical procedure
to transfer-matrix data whose infinite-lattice edge is known analytically (tests/test_edges.py).
"""
from __future__ import annotations

import numpy as np


def attenuation_per_period(Ns, lnT):
    """Least-squares slope of ln T vs N at every frequency.

    ``lnT`` has shape (len(Ns), n_freq). Returns (kappa = -slope, rms residual of the linear fit, in ln units).
    """
    Ns = np.asarray(Ns, float)
    Y = np.asarray(lnT, float)
    A = np.vstack([Ns, np.ones_like(Ns)]).T
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    resid = Y - A @ coef
    return -coef[0], np.sqrt(np.mean(resid ** 2, axis=0))


def edge_from_kappa2(f, kappa, rms, f_seed, side, kappa_window=(0.2, 0.9), max_rms=0.35, min_points=6):
    """Zero crossing of a linear fit of kappa^2 vs f, using points between the seed (mid-gap) and the edge.

    side = 'lo' (edge below the seed) or 'hi'. Points are taken contiguously outward from the seed until kappa leaves the window from above
    (going toward the edge kappa falls); points with kappa below the window or fit rms above ``max_rms`` end the selection.
    Returns (edge, info dict) or (nan, info) if too few points.
    """
    f = np.asarray(f, float)
    s = int(np.argmin(np.abs(f - f_seed)))
    step = -1 if side == "lo" else +1
    idx = []
    i = s
    started = False
    while 0 <= i < len(f):
        k, r = kappa[i], rms[i]
        if not np.isfinite(k) or r > max_rms:
            if started:
                break
        elif k <= kappa_window[1]:
            if k < kappa_window[0]:
                break
            idx.append(i)
            started = True
        i += step
    idx = np.array(idx, int)
    info = dict(n_points=int(len(idx)))
    if len(idx) < min_points:
        return float("nan"), info
    x, y = f[idx], kappa[idx] ** 2
    slope, icpt = np.polyfit(x, y, 1)
    edge = -icpt / slope
    info.update(f_first=float(x.min()), f_last=float(x.max()), slope=float(slope),
                fit_rms=float(np.sqrt(np.mean((y - (slope * x + icpt)) ** 2))))
    return float(edge), info
