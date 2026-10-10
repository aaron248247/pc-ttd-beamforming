"""Phase 6 Part B: system-level error propagation through the (Phase 4) channel transfer functions (pure NumPy).

Perturbation of channel q:   G_q(f) -> G_q(f) a_q exp(i eps_q) exp(+i 2 pi (f - f_ref) dtau_q)
    eps_q   constant phase error (rad)           -> changes ONLY the phase, constant in f
    a_q     amplitude factor 10^(u_q/20)         -> changes ONLY the amplitude
    dtau_q  delay error (in units of a/c)        -> phase slope 2 pi dtau_q in f.  f_ref = 0: uncalibrated (the f0 phase moves too);
                                                    f_ref = f0: calibrated (the intercept at f0 is re-trimmed; pure slope error).
"""
from __future__ import annotations

import numpy as np

from . import array as A


def perturb(G, freqs, eps=None, dtau=None, amp_db=None, f_ref=0.0):
    G = np.asarray(G, complex); f = np.asarray(freqs, float)
    out = G.copy()
    if eps is not None:
        out = out * np.exp(1j * np.asarray(eps, float))[:, None]
    if dtau is not None:
        out = out * np.exp(1j * 2 * np.pi * np.asarray(dtau, float)[:, None] * (f[None, :] - f_ref))
    if amp_db is not None:
        out = out * (10.0 ** (np.asarray(amp_db, float) / 20.0))[:, None]
    return out


def draw(rng, nq, eps_deg=0.0, delay_frac=0.0, amp_db=0.0, tau_rel=None):
    """One independent uniform draw per channel: eps ~ U(-eps_deg, eps_deg) deg, relative delay error e_q ~ U(-delay_frac, delay_frac) times
    tau_rel[q] (the channel's reference delay scale, a/c), amplitude u_q ~ U(-amp_db, amp_db) dB.  Returns (eps rad, dtau a/c, amp_db)."""
    u = lambda w: rng.uniform(-w, w, nq) if w else np.zeros(nq)
    eps = np.deg2rad(u(eps_deg))
    dtau = u(delay_frac) * (np.asarray(tau_rel, float) if tau_rel is not None else 0.0)
    return eps, dtau, u(amp_db)


def percentiles(x):
    x = np.asarray(x, float)
    return dict(median=float(np.median(x)), p90=float(np.percentile(x, 90)), p95=float(np.percentile(x, 95)), worst=float(x.max()))


def monte_carlo(engine, E, G_nom, freqs, theta0, window, p_nom, fixed_max, seedseq, n, level, tau_rel, f_ref=0.0):
    """Monte Carlo at one tolerance point.  level = dict(eps_deg, delay_frac, amp_db).  Returns dict of per-sample arrays:
    max_err, rms_err (deg), gain_loss_db (LOSS, positive = worse: -min over band of 10 log10(p/p_nom), p = peak power incl. any total-amplitude change),
    eff_loss_db (post hoc: the same at fixed input power, i.e. p divided by sum_q |G_q|^2 relative to nominal; isolates taper/phase efficiency from overall level),
    strongest_db (max over band), r_squint."""
    rng = np.random.default_rng(seedseq)
    res = {k: np.empty(n) for k in ("max_err", "rms_err", "gain_loss_db", "eff_loss_db", "strongest_db", "r_squint")}
    S_nom = (np.abs(np.asarray(G_nom)) ** 2).sum(axis=0)
    nq = G_nom.shape[0]
    for i in range(n):
        eps, dtau, adb = draw(rng, nq, level.get("eps_deg", 0.0), level.get("delay_frac", 0.0), level.get("amp_db", 0.0), tau_rel)
        Gp = perturb(G_nom, freqs, eps, dtau, adb, f_ref)
        p = engine.power(Gp, E)
        th, gain, idx = engine.peak(p, theta0, window)
        err = th - theta0
        res["max_err"][i] = np.abs(err).max(); res["rms_err"][i] = np.sqrt(np.mean(err ** 2))
        res["gain_loss_db"][i] = -np.min(10 * np.log10(gain / p_nom))
        res["eff_loss_db"][i] = -np.min(10 * np.log10(gain / p_nom * S_nom / (np.abs(Gp) ** 2).sum(axis=0)))
        res["strongest_db"][i] = engine.strongest_fast(p, idx).max()
        res["r_squint"][i] = res["max_err"][i] / fixed_max
    return res


def budget(levels, p95, limit):
    """Largest half-width (linear interpolation of the P95 curve) with p95 <= limit; None if even the smallest level exceeds it; inf if never exceeded."""
    levels = np.asarray(levels, float); p95 = np.asarray(p95, float)
    if p95[0] > limit:
        return None
    for i in range(1, len(levels)):
        if p95[i] > limit:
            return float(levels[i - 1] + (levels[i] - levels[i - 1]) * (limit - p95[i - 1]) / (p95[i] - p95[i - 1]))
    return float("inf")
