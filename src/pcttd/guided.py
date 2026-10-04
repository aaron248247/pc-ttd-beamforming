"""Guided-mode analysis (pure NumPy): mode projection, 2-port algebra, wave-cascade matrices, delay estimators.

Mode projection (TE: tangential fields Ey, Hz on a plane x = const; Poynting flux S_x = Ey Hz)
    forward mode  phi_f = (Ey_f, Hz_f)          (power flows +x),
    backward mode phi_b = (conj Ey_f, -conj Hz_f)  (time-reversed partner; lossless reciprocal medium).
Lorentz-reciprocity overlap with the total plane field F = (Ey, Hz):
    O_m = integral [ conj(Ey_m) Hz + Ey conj(Hz_m) ] dy,   N_m = same with F = phi_m,   a_m = O_m / N_m.
Because the modes are orthogonal under this form, a_f and a_b are the complex amplitudes of the forward and backward guided mode;
N_b = -N_f, and |a|^2 |N_f| is the modal power.  The profile used is the periodic part u(0, y) at integer x (see w1fdtd), so
a_f(x0) = a_f(0) exp(i k x0): the transfer function between two planes then contains the full propagation phase, H = exp(i k L).

2-port conventions (planes at the device faces): a1 = forward amplitude entering at the left plane, b1 = backward amplitude leaving there,
a2 = backward amplitude entering at the right plane, b2 = forward amplitude leaving there.
    b1 = S11 a1 + S12 a2,   b2 = S21 a1 + S22 a2.
Wave-cascade matrix  [b2; a2] = T [a1; b1]  (T of cascaded 2-ports multiply: T_total = T_B T_A).  A uniform lossless section of length L
with Bloch phase k: T = diag(exp(i k L), exp(-i k L)).
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from . import extract as X


# ----------------------------------------------------------------------------------------------- projection
def resample_periodic(y_src, vals, y_new, H):
    """Cubic periodic interpolation of a profile sampled on [-H/2, H/2) (y_src[0] = -H/2, uniform) onto y_new (complex ok)."""
    n = len(vals)
    dy = H / n
    yy = np.concatenate([y_src, [y_src[0] + H]])
    out = []
    for part in (np.real(vals), np.imag(vals)):
        v = np.concatenate([part, [part[0]]])
        cs = CubicSpline(yy, v, bc_type="periodic")
        yw = (np.asarray(y_new) + H / 2) % H - H / 2
        out.append(cs(yw))
    return out[0] + 1j * out[1]


def project(Ey, Hz, mEy, mHz, w):
    """Forward/backward amplitudes of the plane field (Ey, Hz) on the quadrature ``w`` (weights, e.g. dy)
    onto the forward mode (mEy, mHz). All arrays [nf, ny] (or [ny]).  Returns (a_f, a_b, rel_residual)."""
    Ey, Hz, mEy, mHz = (np.atleast_2d(np.asarray(a, complex)) for a in (Ey, Hz, mEy, mHz))
    w = np.asarray(w, float)
    bEy, bHz = np.conj(mEy), -np.conj(mHz)

    def ov(eE, eH):
        return np.sum(w * (np.conj(eE) * Hz + Ey * np.conj(eH)), axis=-1), np.sum(w * (np.conj(eE) * eH + eE * np.conj(eH)), axis=-1)

    Of, Nf = ov(mEy, mHz)
    Ob, Nb = ov(bEy, bHz)
    af, ab = Of / Nf, Ob / Nb
    rEy = Ey - af[:, None] * mEy - ab[:, None] * bEy
    rHz = Hz - af[:, None] * mHz - ab[:, None] * bHz
    nrm = np.sqrt(np.sum(w * (np.abs(Ey) ** 2 + np.abs(Hz) ** 2), axis=-1))
    res = np.sqrt(np.sum(w * (np.abs(rEy) ** 2 + np.abs(rHz) ** 2), axis=-1)) / np.where(nrm > 0, nrm, np.nan)
    return af, ab, res


# ----------------------------------------------------------------------------------------------- 2-port algebra
def s22_from_unitarity(S11, S21):
    """S22 of a lossless reciprocal 2-port: from S11 conj(S21)... unitarity S11 conj(S12) + S21 conj(S22) = 0 with S12 = S21
    gives S22 = -conj(S11) S21 / conj(S21)."""
    S11, S21 = np.asarray(S11, complex), np.asarray(S21, complex)
    return -np.conj(S11) * S21 / np.conj(S21)


def t_matrix(S11, S21, S22, S12=None):
    """Wave-cascade matrices, shape [nf, 2, 2]: [b2; a2] = T [a1; b1]."""
    S11, S21, S22 = (np.asarray(a, complex) for a in (S11, S21, S22))
    S12 = S21 if S12 is None else np.asarray(S12, complex)
    T = np.empty(S11.shape + (2, 2), complex)
    T[..., 0, 0] = S21 - S11 * S22 / S12
    T[..., 0, 1] = S22 / S12
    T[..., 1, 0] = -S11 / S12
    T[..., 1, 1] = 1.0 / S12
    return T


def s_from_t(T):
    """Inverse of ``t_matrix`` (reciprocal case S12 = S21): returns (S11, S21, S12, S22)."""
    T = np.asarray(T, complex)
    S12 = 1.0 / T[..., 1, 1]
    S11 = -T[..., 1, 0] * S12
    S22 = T[..., 0, 1] * S12
    S21 = T[..., 0, 0] + S11 * S22 / S12
    return S11, S21, S12, S22


def propagate_T(k, L):
    """T of a uniform lossless section: diag(exp(ikL), exp(-ikL)); k, L broadcast."""
    ph = np.exp(1j * np.asarray(k, float) * L)
    T = np.zeros(np.shape(ph) + (2, 2), complex)
    T[..., 0, 0] = ph
    T[..., 1, 1] = 1.0 / ph
    return T


def section_phase_from_lengths(T_long, T_short):
    """Bloch phase theta = k*(N2-N1) (mod 2 pi, sign fixed so that the forward wave is exp(+i theta)) of the section added between two
    otherwise identical devices:  T_long T_short^{-1} = T_R P(dN) T_R^{-1}  has eigenvalues exp(+-i theta), independent of both interfaces.

    The forward eigenvalue is the one whose eigenvector is dominated by the forward component (first entry).
    Returns (theta[nf], eig_modulus[nf]) with theta in (-pi, pi]; unwrap along frequency yourself.
    """
    M = T_long @ np.linalg.inv(T_short)
    th = np.empty(M.shape[0]); mod = np.empty(M.shape[0])
    for i, m in enumerate(M):
        w, v = np.linalg.eig(m)
        j = int(np.argmax(np.abs(v[0]) / np.maximum(np.abs(v[1]), 1e-300)))
        th[i] = np.angle(w[j]); mod[i] = abs(w[j])
    return th, mod


# ----------------------------------------------------------------------------------------------- derivative estimators
def local_poly_derivative(x, y, half: int = 3, order: int = 2):
    """Derivative dy/dx by a least-squares polynomial of ``order`` fitted on 2*half+1 neighboring samples (Savitzky-Golay type).
    End samples use the one-sided window.  Returns an array like y."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); out = np.full(n, np.nan)
    for i in range(n):
        lo = max(0, i - half); hi = min(n, i + half + 1)
        if hi - lo < order + 1:
            continue
        xs = x[lo:hi] - x[i]
        p = np.polyfit(xs, y[lo:hi], order)
        out[i] = p[-2]
    return out


def delay_from_phase(f, phase, method="central", half=3):
    """Group delay +dphi/domega (MEEP_SIGN = +1) from an *unwrapped* phase. method: 'central' (np.gradient) or 'poly' (local quadratic)."""
    f = np.asarray(f, float); w = X.omega(f)
    if method == "central":
        return X.MEEP_SIGN * np.gradient(phase, w)
    return X.MEEP_SIGN * local_poly_derivative(w, phase, half=half)


def unwrap_theta(theta, ref=None):
    """Unwrap a phase in (-pi, pi] along frequency; optional ``ref`` (an approximate unwrapped phase, e.g. from MPB) fixes the 2 pi branch
    at every sample so the result is robust to coarse sampling: theta_unwrapped = theta + 2 pi round((ref - theta)/(2 pi))."""
    theta = np.asarray(theta, float)
    if ref is None:
        return np.unwrap(theta)
    ref = np.asarray(ref, float)
    return theta + 2 * np.pi * np.round((ref - theta) / (2 * np.pi))
