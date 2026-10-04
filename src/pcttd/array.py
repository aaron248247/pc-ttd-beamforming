"""Uniform-linear-array model with frequency-dependent channel transfer functions (Phase 5, pure NumPy).

Conventions (consistent with ``extract.MEEP_SIGN = +1`` and the exp(-i omega t) time dependence of MEEP/MPB)
    * A delay tau multiplies a monochromatic signal by exp(+i omega tau); a propagation length L in a medium gives exp(+i k L).
    * Element n sits at x_n = n d on the array axis; theta is measured from broadside (the normal), positive toward +x.
    * Far field in direction theta:   AF(theta, f) = sum_n H_n(f) exp(-i k0 x_n sin theta),  k0 = 2 pi f / c.
      An element at larger x_n is nearer the observer at positive theta, so the beam steers to theta_0 when
      H_n = exp(+i 2 pi f n d sin(theta_0)/c)  (ideal true-time delay, tau_n = n d sin(theta_0)/c).
    * Isotropic elements (no element pattern); amplitude of H_n carries the channel loss.

Cases
    fixed phase:  H_n = exp(+i 2 pi f0 tau_n)            (frequency-independent phase shifter set at f0)
    ideal TTD:    H_n(f) = exp(+i 2 pi f tau_n)
    measured:     H_n(f) from the FDTD channels; ``trim_phase`` removes one frequency-independent phase per channel (labeled hybrid).
"""
from __future__ import annotations

import numpy as np

C = 1.0


def positions(n, d):
    return np.arange(n) * d


def target_delays(n, d, theta0_deg, c=C):
    return np.arange(n) * d * np.sin(np.deg2rad(theta0_deg)) / c


def fixed_phase_weights(freqs, tau, f0):
    """H[n, f] = exp(+i 2 pi f0 tau_n) for every f."""
    tau = np.asarray(tau, float)
    return np.exp(1j * 2 * np.pi * f0 * tau)[:, None] * np.ones((1, len(freqs)))


def ideal_ttd_weights(freqs, tau):
    return np.exp(1j * 2 * np.pi * np.asarray(tau, float)[:, None] * np.asarray(freqs, float)[None, :])


def trim_phase(H, freqs, tau_target, f0):
    """Phase-intercept trimming (hybrid case): one frequency-independent phase per channel, chosen at f0 only, so that
    arg H_n(f0) - arg H_0(f0) equals 2 pi f0 (tau_target_n - tau_target_0).  Amplitude and frequency dependence are untouched.
    Returns (H_trimmed, delta[n]) with H_trimmed = H exp(-i delta_n)."""
    H = np.asarray(H, complex); freqs = np.asarray(freqs, float)
    i0 = int(np.argmin(np.abs(freqs - f0)))
    ph = np.angle(H[:, i0])
    want = 2 * np.pi * f0 * (np.asarray(tau_target, float) - tau_target[0])
    delta = np.angle(np.exp(1j * ((ph - ph[0]) - want)))
    return H * np.exp(-1j * delta)[:, None], delta


def array_factor(H, freqs, d, sin_theta, c=C):
    """AF[f, s] for channel transfer functions H[n, f] on a grid of sin(theta)."""
    H = np.asarray(H, complex); freqs = np.asarray(freqs, float); st = np.asarray(sin_theta, float)
    n = H.shape[0]
    x = positions(n, d)
    k0 = 2 * np.pi * freqs / c
    ph = np.exp(-1j * k0[:, None, None] * x[None, :, None] * st[None, None, :])       # [f, n, s]
    return np.einsum("nf,fns->fs", H, ph)


def _refine_peak(s, y, i):
    if 0 < i < len(s) - 1:
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            return float(s[i] + 0.5 * (s[1] - s[0]) * (y0 - y2) / den)
    return float(s[i])


def pointing_angle(H, freqs, d, theta0_deg, window=None, n_grid=40001, c=C):
    """Main-lobe pointing angle theta_peak(f) in degrees: the maximum of |AF| over the window |sin(theta) - sin(theta_0)| <= window.

    Default window: half the grating-lobe spacing at the highest frequency, c / (2 f_max d), so that exactly one lobe family is inside.
    Returns (theta_peak[f] in degrees, peak power gain |AF|^2 / N^2 [f], window used)."""
    freqs = np.asarray(freqs, float)
    s0 = np.sin(np.deg2rad(theta0_deg))
    w = c / (2 * freqs.max() * d) if window is None else window
    lo, hi = max(-1.0, s0 - w), min(1.0, s0 + w)
    s = np.linspace(lo, hi, n_grid)
    AF = array_factor(H, freqs, d, s, c)
    p = np.abs(AF) ** 2
    idx = np.argmax(p, axis=1)
    sp = np.array([_refine_peak(s, p[i], idx[i]) for i in range(len(freqs))])
    return np.rad2deg(np.arcsin(np.clip(sp, -1, 1))), p[np.arange(len(freqs)), idx] / H.shape[0] ** 2, w


def pattern(H_f, f, d, n_grid=4001, c=C):
    """|AF|^2 / N^2 over the full visible region for one frequency; returns (theta_deg, power)."""
    th = np.linspace(-90, 90, n_grid)
    AF = array_factor(np.asarray(H_f, complex)[:, None], [f], d, np.sin(np.deg2rad(th)), c)[0]
    return th, np.abs(AF) ** 2 / len(H_f) ** 2


def half_power_beamwidth(theta_deg, p, i_peak):
    """Full width (degrees) between the -3 dB points around the peak index; nan if the pattern does not fall below -3 dB on a side."""
    half = p[i_peak] / 2.0
    j = i_peak
    while j > 0 and p[j] > half: j -= 1
    k = i_peak
    while k < len(p) - 1 and p[k] > half: k += 1
    if j == 0 or k == len(p) - 1:
        return float("nan")
    return float(theta_deg[k] - theta_deg[j])


def r_squint(theta_phc, theta_phase, theta0_deg, mask=None):
    """R_squint = max_f |theta_PhC - theta_0| / max_f |theta_phase - theta_0| over the (boolean) frequency mask."""
    a = np.abs(np.asarray(theta_phc, float) - theta0_deg)
    b = np.abs(np.asarray(theta_phase, float) - theta0_deg)
    if mask is not None:
        a, b = a[mask], b[mask]
    return float(a.max() / b.max())
