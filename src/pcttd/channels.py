"""Phase 4: four-channel true-time-delay family (pure NumPy).

Conventions (a = c = 1, f = a/lambda, omega = 2 pi f, MEEP exp(-i omega t), MEEP_SIGN = +1): a delay tau multiplies a monochromatic signal by
exp(+i omega tau), tau_g = +d arg H / d omega.  Channel m is a directly simulated finite device of N_m = N0 + m dN periods; H_m(f) = S21_m is the
complete transfer function between the common reference planes (it contains the whole propagation phase; nothing is de-embedded).
The ideal channel family relative to channel 0 is  H_m / H_0 = exp(+i 2 pi f m dtau_target).

Raw channels use H_m as simulated.  The phase-intercept-trimmed channels multiply each H_m by exp(-i delta_m), one constant per channel chosen at f0
only (``trim_constants``); this changes neither |H_m| nor the frequency dependence of the phase, hence not the group delay.
"""
from __future__ import annotations

import numpy as np

from . import guided as GD
from . import ttdmetrics as TM


def channel_lengths(N0, dN, n=4):
    """N_m = N0 + m dN, m = 0..n-1."""
    if N0 < 1 or dN < 1:
        raise ValueError("N0 and dN must be positive integers")
    return [int(N0 + m * dN) for m in range(n)]


def ideal_channels(f, dtau_target, n=4, amp=1.0):
    """Ideal TTD family H_m(f) = amp * exp(+i 2 pi f m dtau_target), shape [n, nf]."""
    f = np.asarray(f, float)
    m = np.arange(n)[:, None]
    return amp * np.exp(1j * 2 * np.pi * f[None, :] * m * dtau_target)


def relative_phase(H, ref_per_period=None, dN=1):
    """Unwrapped arg(H_m / H_0), shape [n, nf].  The 2 pi branch of every sample is fixed by ``ref_per_period`` (an approximate unwrapped
    phase advance per period, e.g. MPB -2 pi kappa; reference for channel m is m dN times it); without it the phase is unwrapped along f."""
    H = np.asarray(H, complex)
    out = np.empty(H.shape)
    for m in range(H.shape[0]):
        w = np.angle(H[m] / H[0])
        if ref_per_period is None:
            out[m] = np.unwrap(w)
        else:
            out[m] = GD.unwrap_theta(w, ref=np.asarray(ref_per_period, float) * m * dN)
    return out


def delays(f, phi, method="central", half=3):
    """tau_g = +d phi / d omega for each row of the unwrapped phase ``phi`` [n, nf]."""
    return np.array([GD.delay_from_phase(f, p, method=method, half=half) for p in np.atleast_2d(phi)])


def incremental_delays(f, H, ref_per_period=None, dN=1, method="central", half=3):
    """Adjacent increments dtau_{m+1,m}(f) = d arg(H_{m+1}/H_m)/d omega, shape [n-1, nf] (branch fixed per pair by ``ref_per_period`` * dN)."""
    H = np.asarray(H, complex)
    rows = []
    for m in range(H.shape[0] - 1):
        w = np.angle(H[m + 1] / H[m])
        p = np.unwrap(w) if ref_per_period is None else GD.unwrap_theta(w, ref=np.asarray(ref_per_period, float) * dN)
        rows.append(GD.delay_from_phase(f, p, method=method, half=half))
    return np.array(rows)


def delay_error(dtau_m, dtau_target):
    """Relative delays Dtau_m(f) [n, nf] (row 0 is zero) against m * dtau_target: returns (absolute, normalized) with normalized NaN for m = 0."""
    dtau_m = np.asarray(dtau_m, float)
    m = np.arange(dtau_m.shape[0])[:, None]
    err = dtau_m - m * dtau_target
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(m > 0, err / (m * dtau_target), np.nan)
    return err, rel


def raw_phase_error(f, phi_rel, dtau_target, f0):
    """epsilon_phi_m(f) = phi_rel_m(f) - 2 pi f m dtau_target, placed on the 2 pi branch with epsilon_phi_m(f0) in (-pi, pi].  Shape [n, nf]."""
    f = np.asarray(f, float)
    phi_rel = np.asarray(phi_rel, float)
    m = np.arange(phi_rel.shape[0])[:, None]
    e = phi_rel - 2 * np.pi * f[None, :] * m * dtau_target
    i0 = int(np.argmin(np.abs(f - f0)))
    shift = np.angle(np.exp(1j * e[:, i0])) - e[:, i0]          # multiple of 2 pi taking e(f0) into (-pi, pi]
    return e + (2 * np.pi * np.round(shift / (2 * np.pi)))[:, None]


def trim_constants(eps_raw, f, f0):
    """delta_m = raw phase error at f0 (frequency-independent trim per channel); delta_0 = 0 by construction."""
    i0 = int(np.argmin(np.abs(np.asarray(f, float) - f0)))
    return np.asarray(eps_raw, float)[:, i0].copy()


def apply_trim(H, delta):
    """H_m exp(-i delta_m): one constant phase per channel; amplitude and frequency dependence untouched."""
    return np.asarray(H, complex) * np.exp(-1j * np.asarray(delta, float))[:, None]


def amplitude_imbalance_db(H):
    """20 log10(max_m |H_m| / min_m |H_m|) per frequency (AMPLITUDE ratio in dB, not a power ratio)."""
    a = np.abs(np.asarray(H, complex))
    return 20 * np.log10(a.max(axis=0) / a.min(axis=0))


def common_usable_mask(masks):
    """Logical AND of per-channel / per-criterion boolean masks."""
    out = np.ones_like(np.asarray(masks[0], bool))
    for m in masks:
        out &= np.asarray(m, bool)
    return out


def widest_run(mask, f):
    """(i_low, i_high) of the widest contiguous run of True in ``mask`` (width measured in f); None if there is none."""
    runs = TM.contiguous_runs(mask)
    if not runs:
        return None
    f = np.asarray(f, float)
    return max(runs, key=lambda r: f[r[1]] - f[r[0]])


def cascade_predict(S11a, S21a, S11b, S21b, steps):
    """Predict the 2-port of the device ``steps`` periods longer than device b (steps * (Nb - Na) added periods in total, with Nb - Na the length difference
    of a and b): T_c = [T_b T_a^-1]^steps T_b.  Returns (S11, S21) of the prediction.  ``steps`` is a positive integer."""
    Ta = GD.t_matrix(S11a, S21a, GD.s22_from_unitarity(S11a, S21a))
    Tb = GD.t_matrix(S11b, S21b, GD.s22_from_unitarity(S11b, S21b))
    M = Tb @ np.linalg.inv(Ta)
    T = Tb
    for _ in range(int(steps)):
        T = M @ T
    S11, S21, _, _ = GD.s_from_t(T)
    return S11, S21
