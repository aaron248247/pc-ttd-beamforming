"""Phase 6 hierarchical (subarray-level TTD) array model on top of ``pcttd.array`` / ``pcttd.beam`` (pure NumPy).

Element (q, p) of subarray q (q = 0..nq-1, p = 0..msub-1) sits at x = (q msub + p) d_e and carries the weight
    W_(q,p)(f) = G_q(f) * U_p(f),     G_q = subarray channel transfer function (PhC, ideal TTD, fixed phase),
                                       U_p = intra-subarray weight (normally the fixed phase exp(+i 2 pi f0 p d_e sin(theta0))).
The far field factorizes exactly:  AF(s, f) = [sum_q G_q exp(-i k0 q msub d_e s)] [sum_p U_p exp(-i k0 p d_e s)]  (s = sin(theta)),
so each subarray is NOT treated as a point source; the element factor E(s, f) is the full intra-subarray sum.
Power is |AF|^2 / N_el^2, N_el = nq msub, as in ``pcttd.array``.  Conventions (exp(-i omega t); delay tau -> exp(+i omega tau)) are those of ``pcttd.array``.
"""
from __future__ import annotations

import numpy as np

from . import array as A
from . import beam as B

C = A.C


def tau_required(msub, d_e, theta0_deg, c=C):
    """Ideal delay between adjacent subarray centers/origins: msub d_e sin(theta0) / c."""
    return msub * d_e * np.sin(np.deg2rad(theta0_deg)) / c


def mismatch(tau_phc, tau_req, f0):
    """(rho, rho - 1, f0 (tau_phc - tau_req) in carrier cycles)."""
    return tau_phc / tau_req, tau_phc / tau_req - 1.0, f0 * (tau_phc - tau_req)


def intra_fixed(msub, d_e, theta0_deg, f0, nf, c=C):
    """U[p, f] = exp(+i 2 pi f0 p d_e sin(theta0)/c), frequency independent."""
    u = np.exp(1j * 2 * np.pi * f0 * np.arange(msub) * d_e * np.sin(np.deg2rad(theta0_deg)) / c)
    return np.repeat(u[:, None], nf, axis=1)


def intra_ttd(msub, d_e, theta0_deg, freqs, c=C):
    """U[p, f] = exp(+i 2 pi f p d_e sin(theta0)/c) (ideal TTD inside the subarray)."""
    tau = np.arange(msub) * d_e * np.sin(np.deg2rad(theta0_deg)) / c
    return np.exp(1j * 2 * np.pi * tau[:, None] * np.asarray(freqs, float)[None, :])


def ideal_channels(freqs, nq, tau_req):
    """G_q = exp(+i 2 pi f q tau_req): exact TTD between subarrays."""
    return A.ideal_ttd_weights(freqs, np.arange(nq) * tau_req)


def fixed_channels(freqs, nq, tau_req, f0):
    return A.fixed_phase_weights(freqs, np.arange(nq) * tau_req, f0)


def retrim(G, freqs, tau_req, f0):
    """Phase-intercept re-trim to the REQUIRED step: one constant per channel so arg G_q(f0) - arg G_0(f0) = 2 pi f0 q tau_req.
    Group delay and amplitude untouched.  Returns (G_trim, delta)."""
    return A.trim_phase(G, freqs, np.arange(G.shape[0]) * tau_req, f0)


def hier_weights(G, U):
    """Element-level weights W[(q msub + p), f] = G[q, f] U[p, f]."""
    G = np.asarray(G, complex); U = np.asarray(U, complex)
    return (G[:, None, :] * U[None, :, :]).reshape(G.shape[0] * U.shape[0], G.shape[1])


def window_halfwidth(fmax, msub, d_e, c=C):
    """Design-lobe window in sin(theta): half the subarray-level grating-lobe spacing at fmax."""
    return c / (2 * fmax * msub * d_e)


class Hier:
    """Fast exact evaluator of the full far-field power on a uniform sin(theta) grid for the factorized hierarchical array."""

    def __init__(self, freqs, nq, msub, d_e, n=40001, c=C):
        self.f = np.asarray(freqs, float); self.nq, self.msub, self.d_e, self.c = nq, msub, d_e, c
        self.s = np.linspace(-1.0, 1.0, n)
        self.k = 2 * np.pi * self.f / c
        self.Pq = np.exp(-1j * self.k[:, None, None] * (np.arange(nq) * msub * d_e)[None, :, None] * self.s[None, None, :])   # [f, q, s]
        self.N = nq * msub

    def element_factor(self, U):
        """E[f, s] = sum_p U[p, f] exp(-i k p d_e s)."""
        U = np.asarray(U, complex)
        E = np.empty((len(self.f), len(self.s)), complex)
        p = np.arange(self.msub)
        for i, k in enumerate(self.k):
            E[i] = U[:, i] @ np.exp(-1j * k * p[:, None] * self.d_e * self.s[None, :])
        return E

    def af(self, G, E):
        return np.einsum("qf,fqs->fs", np.asarray(G, complex), self.Pq) * E

    def power(self, G, E):
        return np.abs(self.af(G, E)) ** 2 / self.N ** 2

    # ---- metrics ---------------------------------------------------------------------------------------------------------------------------
    def _local_max(self, p):
        m = np.zeros(p.shape, bool)
        m[:, 1:-1] = (p[:, 1:-1] > p[:, :-2]) & (p[:, 1:-1] >= p[:, 2:])
        m[:, 0] = p[:, 0] > p[:, 1]
        m[:, -1] = p[:, -1] > p[:, -2]
        return m

    def peak(self, p, theta0_deg, window):
        """Design-lobe peak inside |s - s0| <= window: (theta_peak deg, peak power, index array)."""
        s0 = np.sin(np.deg2rad(theta0_deg))
        lo, hi = np.searchsorted(self.s, max(-1.0, s0 - window)), np.searchsorted(self.s, min(1.0, s0 + window), side="right")
        idx = lo + np.argmax(p[:, lo:hi], axis=1)
        r = np.arange(p.shape[0])
        i = np.clip(idx, 1, p.shape[1] - 2)
        y0, y1, y2 = p[r, i - 1], p[r, i], p[r, i + 1]
        den = y0 - 2 * y1 + y2
        ds = self.s[1] - self.s[0]
        sp = np.where((idx == i) & (den != 0), self.s[i] + 0.5 * ds * (y0 - y2) / np.where(den == 0, 1, den), self.s[idx])
        return np.rad2deg(np.arcsin(np.clip(sp, -1, 1))), p[r, idx], idx

    def lobes(self, p, idx):
        """Strongest local maximum other than the design lobe (the local maximum nearest the peak index), dB below the design peak, and its angle."""
        lm = self._local_max(p)
        out_db = np.empty(len(p)); out_th = np.empty(len(p)); n3 = np.empty(len(p), int); kd = np.empty(len(p), int)
        for i in range(len(p)):
            m = np.flatnonzero(lm[i])
            k = m[np.argmin(np.abs(m - idx[i]))]
            kd[i] = k
            pk = p[i, k]
            o = m[m != k]
            if len(o):
                j = o[np.argmax(p[i, o])]
                out_db[i] = 10 * np.log10(max(p[i, j], 1e-300) / pk); out_th[i] = np.rad2deg(np.arcsin(self.s[j]))
            else:
                out_db[i], out_th[i] = -np.inf, np.nan
            mi = m[(m > 0) & (m < p.shape[1] - 1)]
            n3[i] = int(np.sum(p[i, mi] >= pk / 2))
        return out_db, out_th, n3, kd

    def strongest_fast(self, p, idx):
        """Vectorized strongest-other-lobe level (dB re the design peak) for Monte Carlo; the design lobe is the local maximum at ``idx``."""
        lm = self._local_max(p)
        r = np.arange(len(p))
        pk = p[r, idx]
        lm[r, idx] = False
        mx = np.where(lm, p, 0.0).max(axis=1)
        return 10 * np.log10(np.maximum(mx, 1e-300) / pk)

    def replica_levels(self, p, theta_peak, half=0.03):
        """Level (dB re the design peak) of the strongest power within +-half in sin(theta) of each analytic subarray-level replica of the design lobe,
        sin = sin(theta_peak) + q c / (f nsub_pitch), q != 0, |sin| <= 1 (pitch = msub d_e).  Returns (level_db[f] (-inf if no visible replica), n_visible[f])."""
        lev = np.full(len(p), -np.inf); nvis = np.zeros(len(p), int)
        pitch = self.msub * self.d_e
        for i in range(len(p)):
            sp = np.sin(np.deg2rad(theta_peak[i]))
            pk = p[i, np.argmin(np.abs(self.s - sp))]
            q = np.arange(-16, 17); q = q[q != 0]
            sg = sp + q * self.c / (self.f[i] * pitch)
            sg = sg[np.abs(sg) <= 1 + 1e-9]
            nvis[i] = len(sg)
            for s_ in sg:
                m = np.abs(self.s - s_) <= half
                lev[i] = max(lev[i], 10 * np.log10(max(p[i, m].max(), 1e-300) / pk))
        return lev, nvis

    def analyze(self, G, U, theta0_deg, window):
        """Per-frequency design-lobe metrics, same keys as ``beam.analyze_case`` (plus the angle of the strongest other lobe)."""
        E = self.element_factor(U)
        p = self.power(G, E)
        th, gain, idx = self.peak(p, theta0_deg, window)
        db, tho, n3, kd = self.lobes(p, idx)
        hp = np.array([B._hpbw(self.s, p[i], kd[i]) for i in range(len(p))])
        edge = np.stack([10 * np.log10(np.maximum(p[:, 0], 1e-300) / p[np.arange(len(p)), kd]),
                         10 * np.log10(np.maximum(p[:, -1], 1e-300) / p[np.arange(len(p)), kd])], axis=1)
        gl, nvis = self.replica_levels(p, th)
        return dict(theta_peak=th, gain=gain, hpbw=hp, strongest_db=db, strongest_theta=tho, n_3db=n3, edge_db=edge, grating_db=gl, n_replicas=nvis)

    def pattern(self, G_f, U_f, i):
        """Full normalized power over the grid at frequency index i for subarray weights G_f[q], intra weights U_f[p]; returns (theta_deg, power)."""
        p_ = np.arange(self.msub)
        E = U_f @ np.exp(-1j * self.k[i] * p_[:, None] * self.d_e * self.s[None, :])
        AF = np.asarray(G_f, complex) @ self.Pq[i] * E
        return np.rad2deg(np.arcsin(self.s)), np.abs(AF) ** 2 / self.N ** 2
