"""Phase and group-delay extraction from complex transfer functions.

Units: a = 1, c = 1; f = omega a/(2 pi c) is the frequency in cycles, omega = 2 pi f; delays are in a/c.

Fourier-sign convention (derived here, then verified numerically against FDTD in Phase 1A)
----------------------------------------------------------------------------------------
MEEP time-steps fields with time dependence exp(-i omega t) and accumulates the DFT
    F(omega) = integral f(t) exp(+i omega t) dt ,
which picks out the exp(-i omega t) component of a real signal. A right-going wave is
    E(x, t) = Re[ E0 exp(i (k x - omega t)) ]  ->  F(omega, x) = E0 exp(i k x).
For the transfer function between two reference planes separated by length L,
    H(omega) = F(omega, x + L) / F(omega, x) = exp(i k L),        phi = arg H = k L.
The group delay of a wavepacket is tau_g = L / v_g = L dk/d omega = + d phi / d omega.
With the opposite (exp(+i omega t)) convention the sign flips: tau_g = - d phi / d omega.
The ``sign`` argument below selects the convention; the project sign is ``MEEP_SIGN = +1``, and that is
asserted against vacuum and homogeneous-dielectric FDTD runs with known positive delay (tests/test_fdtd_1d.py),
not taken on faith.
"""
from __future__ import annotations

import warnings

import numpy as np

MEEP_SIGN = +1          # tau_g = MEEP_SIGN * d(arg H)/d(omega); verified numerically in Phase 1A
TWO_PI = 2.0 * np.pi


def omega(f):
    return TWO_PI * np.asarray(f, float)


def max_unambiguous_delay(df: float) -> float:
    """Largest |delay| whose phase ramp changes by less than pi between adjacent frequency samples.

    A delay tau adds phase 2 pi tau df per sample; unwrapping is unambiguous only if that is < pi,
    i.e. tau < 1 / (2 df).
    """
    return 1.0 / (2.0 * df)


def unwrap_checked(phase_wrapped, max_step: float = 0.9 * np.pi):
    """np.unwrap plus a diagnostic.

    Returns (unwrapped, ok, max_wrapped_step). ``ok`` is False when any wrapped sample-to-sample phase step
    exceeds ``max_step``. This is a *necessary* sanity check only: a delay that aliases to a small apparent step passes
    it, so use ``assert_sampling_adequate`` with the known maximum delay to guarantee a valid unwrap.
    Samples that are NaN (masked) are skipped, and the unwrapping restarts across them.
    """
    p = np.asarray(phase_wrapped, float)
    out = np.full_like(p, np.nan)
    good = np.isfinite(p)
    if not good.any():
        return out, True, 0.0
    idx = np.nonzero(good)[0]
    out[idx] = np.unwrap(p[idx])
    steps = np.abs(np.angle(np.exp(1j * np.diff(p[idx]))))
    mx = float(steps.max()) if steps.size else 0.0
    return out, bool(mx <= max_step), mx


def group_delay_phase(f, H, sign: int = MEEP_SIGN, mask_below: float | None = None):
    """tau_g = sign * d arg(H) / d omega from the unwrapped phase (second-order central differences).

    ``mask_below``: |H| values below this are set to NaN before unwrapping (phase is meaningless there).
    Returns (tau, ok) where ok is the unwrapping-validity flag.
    """
    f = np.asarray(f, float)
    H = np.asarray(H, complex)
    ph = np.angle(H)
    if mask_below is not None:
        ph = np.where(np.abs(H) >= mask_below, ph, np.nan)
    un, ok, _ = unwrap_checked(ph)
    tau = np.full_like(f, np.nan)
    g = np.isfinite(un)
    if g.sum() >= 3:
        tau[g] = sign * np.gradient(un[g], omega(f)[g])
    return tau, ok


def group_delay_ratio(f, H, sign: int = MEEP_SIGN):
    """Unwrapping-free group delay: tau = sign * arg(H[i+1] conj(H[i-1])) / (omega[i+1] - omega[i-1]).

    The principal value of the phase of the neighbor ratio equals the phase *difference* across the two-sample stencil
    as long as that difference is below pi (i.e. 2 * d_omega * |tau| < pi), so no global unwrapping is needed. It is exact
    for phase that is linear or quadratic in omega. End points use one-sided differences.
    """
    f = np.asarray(f, float)
    H = np.asarray(H, complex)
    w = omega(f)
    out = np.full_like(f, np.nan)
    ratio = H[2:] * np.conj(H[:-2])
    out[1:-1] = sign * np.angle(ratio) / (w[2:] - w[:-2])
    out[0] = sign * np.angle(H[1] * np.conj(H[0])) / (w[1] - w[0])
    out[-1] = sign * np.angle(H[-1] * np.conj(H[-2])) / (w[-1] - w[-2])
    return out


def group_delay_derivative(f, H, sign: int = MEEP_SIGN):
    """tau = sign * Im[(dH/d omega)/H] with a central-difference dH/d omega.

    Needs no unwrapping, but differentiating the oscillating H has truncation error ~ tau * (d_omega * tau)^2 / 6,
    so it is only accurate when d_omega * tau << 1 (much stricter than the unwrapping limit). Kept as a cross-check
    whose error scaling is tested; ``group_delay_ratio`` is the unwrapping-free method used for results.
    """
    f = np.asarray(f, float)
    H = np.asarray(H, complex)
    dH = np.gradient(H, omega(f))
    with np.errstate(divide="ignore", invalid="ignore"):
        return sign * np.imag(dH / H)


def assert_sampling_adequate(f, tau_max: float, safety: float = 0.8, span: int = 1):
    """Raise unless the frequency grid is fine enough to unwrap a phase whose group delay may reach ``tau_max``.

    Aliasing cannot be detected from the wrapped phase alone (a large delay can alias to a small apparent step, e.g.
    a step of exactly 4 pi looks like zero), so the requirement is imposed a priori from the physical size of the
    structure: max(d_f) <= safety / (2 tau_max * span). ``span`` is the number of sample steps the method's stencil
    covers (1 for np.unwrap-based extraction, 2 for ``group_delay_ratio``).
    """
    df = float(np.max(np.diff(np.asarray(f, float))))
    limit = safety * max_unambiguous_delay(df) / span
    if tau_max > limit:
        raise ValueError(f"frequency spacing {df:g} cannot unwrap delays up to {tau_max:g} "
                         f"(limit {limit:g}); use spacing <= {safety / (2 * tau_max * span):g}")
    return True


def band_mask(f, fmin, fmax):
    f = np.asarray(f, float)
    return (f >= fmin) & (f <= fmax)


def fit_delay_vs_length(lengths, delays):
    """Least-squares line delay = slope * L + intercept. Returns (slope, intercept, max_abs_residual)."""
    L = np.asarray(lengths, float)
    d = np.asarray(delays, float)
    A = np.vstack([L, np.ones_like(L)]).T
    (slope, icpt), *_ = np.linalg.lstsq(A, d, rcond=None)
    return float(slope), float(icpt), float(np.max(np.abs(A @ np.array([slope, icpt]) - d)))


def observed_order(h, err):
    """Least-squares slope of log(err) vs log(h): the measured convergence order (not an assumed one).

    ``h`` grid spacings (1/resolution); ``err`` positive error magnitudes. Returns (order, intercept_log).
    """
    h = np.asarray(h, float)
    e = np.asarray(err, float)
    ok = np.isfinite(e) & (e > 0)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    p, c = np.polyfit(np.log(h[ok]), np.log(e[ok]), 1)
    return float(p), float(c)


# ----------------------------------------------------------------------------- Yee-grid prediction
def yee_1d_delay_phase(f, L, n, dx, dt):
    """Exact numerical-dispersion prediction for a homogeneous 1-D Yee grid (c = 1).

    Discrete dispersion: sin(k dx / 2) = (n dx/dt) sin(omega dt / 2).  Returns (phase kL, group delay L dk/d omega).
    This is what FDTD *should* give in a homogeneous medium: the residual of a measurement against it
    tests the extraction, while the deviation of the prediction from the ideal n L is the (known) numerical dispersion.
    """
    w = omega(f)
    s = n * (dx / dt) * np.sin(w * dt / 2.0)
    k = (2.0 / dx) * np.arcsin(s)
    dk = (2.0 / dx) * (n * (dx / dt) * (dt / 2.0) * np.cos(w * dt / 2.0)) / np.sqrt(1.0 - s * s)
    return k * L, L * dk
