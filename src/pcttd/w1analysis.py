"""Phase 2+ analysis of finite-W1 runs: modal S-parameters from saved plane fields, multi-length de-embedding, delay estimators.

Pure NumPy except ``modal_sparams`` which needs the MPB mode profiles (``w1modes``) and is passed them already evaluated.

Conventions: a = c = 1, f = a/lambda, omega = 2 pi f, MEEP_SIGN = +1 (tau_g = +d arg H / d omega), the forward guided wave has
Bloch wavenumber k_fwd = -kappa with kappa the MPB wavenumber of the decreasing Hz-even band; the *phase advance per period* used here is
theta(f) = k_fwd-convention phase of the section transfer, i.e. S21 ~ exp(+i theta N), d theta/d omega = n_g (a/c).
"""
from __future__ import annotations

import numpy as np

from . import extract as X
from . import guided as GD


def modal_sparams(raw: dict, mode_y, mode_Ey, mode_Hz, H):
    """S21, S11 and projection residuals from a ``w1fdtd.run_w1`` result.

    ``mode_*`` are the forward-mode profiles on ``mode_y`` (one per frequency, shape [nf, ny]) over the supercell of height ``H``.
    Returns dict(S21, S11, res_in, res_out, ab_out) with ab_out = backward amplitude at the output plane divided by af_in (spurious return)."""
    out = {}
    amps = {}
    for name in ("in", "out"):
        y = raw[name + "_y"]
        mEy = np.array([GD.resample_periodic(mode_y, mode_Ey[i], y, H) for i in range(len(mode_Ey))])
        mHz = np.array([GD.resample_periodic(mode_y, mode_Hz[i], y, H) for i in range(len(mode_Hz))])
        amps[name] = GD.project(raw[name + "_Ey"], raw[name + "_Hz"], mEy, mHz, raw[name + "_w"])
    af_in, ab_in, r_in = amps["in"]
    af_out, ab_out, r_out = amps["out"]
    return dict(S21=af_out / af_in, S11=ab_in / af_in, res_in=r_in, res_out=r_out, ab_out=ab_out / af_in)


def unitarity_defect(S11, S21):
    """1 - |S11|^2 - |S21|^2 (zero for a lossless reciprocal two-port whose faces lie in the single-mode lead)."""
    return 1.0 - np.abs(S11) ** 2 - np.abs(S21) ** 2


def cascade_theta(S11_by_N: dict, S21_by_N: dict, N_a: int, N_b: int):
    """Per-period phase advance from two lengths by the exact wave-cascade method (independent of both interfaces).
    Returns (theta_total[nf] in (-pi, pi], modulus[nf]) for the added length N_b - N_a."""
    Ts = {}
    for N in (N_a, N_b):
        S11, S21 = S11_by_N[N], S21_by_N[N]
        Ts[N] = GD.t_matrix(S11, S21, GD.s22_from_unitarity(S11, S21))
    return GD.section_phase_from_lengths(Ts[N_b], Ts[N_a])


def unwrap_to_reference(theta, theta_ref):
    """Branch-fix a wrapped phase to the nearest 2 pi branch of a reference (e.g. MPB) phase."""
    return GD.unwrap_theta(theta, ref=theta_ref)


def linear_length_fit(N, phases, weights=None):
    """Fit phase_i(f) = c(f) + s(f) N_i for every frequency. ``phases`` [nN, nf] (already unwrapped consistently).
    Returns (slope[nf], intercept[nf], rms_residual[nf]); with 3+ lengths the residual tests linear length scaling."""
    N = np.asarray(N, float); P = np.asarray(phases, float)
    A = np.vstack([N, np.ones_like(N)]).T
    if weights is not None:
        A = A * np.asarray(weights)[:, None]; P = P * np.asarray(weights)[:, None]
    sol, *_ = np.linalg.lstsq(A, P, rcond=None)
    res = P - A @ sol
    dof = max(len(N) - 2, 1)
    return sol[0], sol[1], np.sqrt(np.sum(res ** 2, axis=0) / dof)


def ng_from_theta(f, theta_per_period, method="central", half=3):
    """n_g = d theta_per_period / d omega (a = c = 1), MEEP sign convention. method as in ``guided.delay_from_phase``."""
    return GD.delay_from_phase(f, theta_per_period, method=method, half=half)


def subsample(f, arrays, step):
    """Take every ``step``-th frequency of the fine grid (used to test the frequency-sampling sensitivity of the derivative)."""
    return np.asarray(f)[::step], [np.asarray(a)[::step] for a in arrays]


def ng_vs_reference(f, ng, ng_ref, mask=None):
    """Absolute and relative discrepancy, with summary statistics over ``mask``."""
    d = np.asarray(ng) - np.asarray(ng_ref)
    rel = d / np.asarray(ng_ref)
    m = np.ones_like(rel, bool) if mask is None else mask
    return dict(abs=d, rel=rel, max_abs_rel=float(np.nanmax(np.abs(rel[m]))), rms_rel=float(np.sqrt(np.nanmean(rel[m] ** 2))),
                median_rel=float(np.nanmedian(rel[m])))


def pulse_centroid_delay(t, e1, e2):
    """Centroid delay between two time traces: first moment of |e|^2 (energy centroid) difference."""
    t = np.asarray(t, float)
    w1, w2 = np.abs(e1) ** 2, np.abs(e2) ** 2
    return float(np.sum(t * w2) / np.sum(w2) - np.sum(t * w1) / np.sum(w1))
