"""Pure-Python transfer-matrix method for a lossless 1-D dielectric multilayer at normal incidence.

Convention: time dependence exp(-i omega t) (same as MEEP). In medium j (index n_j, real)
    E(x) = a_j exp(+i n_j k0 (x - x_j)) + b_j exp(-i n_j k0 (x - x_j)),   k0 = 2 pi f,   H_y ~ n_j (a_j e^{..} - b_j e^{..}),
with x_j the left face of the layer, so a is right-going and b left-going. Continuity of E and H at the
interface between media j and j+1 (rho = n_{j+1}/n_j) gives, for amplitudes referenced at that face,
    [a;b]_j = D [a;b]_{j+1},   D = (1/2) [[1+rho, 1-rho],[1-rho, 1+rho]],
and crossing a layer of thickness d moves the reference plane: [a;b]_left = diag(e^{-i delta}, e^{+i delta}) [a;b]_right,
delta = n k0 d. Chaining and requiring no incident wave from the right (b_s = 0):
    [a0; b0] = M [a_s; 0]   ->   t = a_s/a0 = 1/M00,   r = b0/a0 = M10/M00.
t is referenced face-to-face (phase measured at the entrance face for the input, at the exit face for the output);
r is referenced to the entrance face. For a vacuum "layer" of thickness d, t = exp(i k0 d), so arg t = +k d and
tau_g = +d arg(t)/d omega = d/c, matching the MEEP convention derived in extract.py.

Power: T = (n_s/n_0) |t|^2, R = |r|^2, and T + R = 1 for lossless layers.
"""
from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def _chain(f, layers, n_in, n_out):
    f = np.atleast_1d(np.asarray(f, float))
    M = np.zeros((f.size, 2, 2), complex)
    M[:] = np.eye(2)
    ns = [n_in] + [n for n, _ in layers] + [n_out]
    ds = [0.0] + [d for _, d in layers]
    for j in range(len(ns) - 1):
        rho = ns[j + 1] / ns[j]
        D = 0.5 * np.array([[1 + rho, 1 - rho], [1 - rho, 1 + rho]], complex)
        delta = ns[j] * TWO_PI * f * ds[j]
        P = np.zeros((f.size, 2, 2), complex)
        P[:, 0, 0] = np.exp(-1j * delta)
        P[:, 1, 1] = np.exp(+1j * delta)
        # reference plane of layer j is its left face; propagate to its right face then cross interface j|j+1
        M = M @ (P @ D)
    return M


def transmission_reflection(f, layers, n_in=1.0, n_out=1.0):
    """Complex t (face-to-face), complex r (entrance face), power T and R. ``layers``: list of (n, thickness)."""
    M = _chain(f, layers, n_in, n_out)
    t = 1.0 / M[:, 0, 0]
    r = M[:, 1, 0] / M[:, 0, 0]
    T = (n_out / n_in) * np.abs(t) ** 2
    R = np.abs(r) ** 2
    return t, r, T, R


def total_thickness(layers):
    return float(sum(d for _, d in layers))


def vacuum_referenced_t(f, layers, n_in=1.0, n_out=1.0):
    """t multiplied by exp(-i k0 D): the transmission relative to propagating through the same total thickness D of
    vacuum, which is what an FDTD normalization run (same monitor planes, structure removed) measures when n_in=n_out=1."""
    f = np.atleast_1d(np.asarray(f, float))
    t, *_ = transmission_reflection(f, layers, n_in, n_out)
    return t * np.exp(-1j * TWO_PI * f * total_thickness(layers))


def tmm_group_delay(f, layers, n_in=1.0, n_out=1.0, referenced=True, h=1e-6):
    """Group delay d arg(t)/d omega of the transmission (excess over vacuum if ``referenced``) by central difference
    of the complex t on +-h (frequency units); analytic phase, no unwrapping needed."""
    f = np.atleast_1d(np.asarray(f, float))
    fn = vacuum_referenced_t if referenced else (lambda ff, L, a, b: transmission_reflection(ff, L, a, b)[0])
    tp, tm, t0 = fn(f + h, layers, n_in, n_out), fn(f - h, layers, n_in, n_out), fn(f, layers, n_in, n_out)
    dt = (tp - tm) / (TWO_PI * 2 * h)
    return np.imag(dt / t0)


# ----------------------------------------------------------------------------- periodic structures
def bilayer(n_a, d_a, n_b, d_b, N):
    return [(n_a, d_a), (n_b, d_b)] * int(N)


def quarter_wave_bilayer(n_lo, n_hi, period=1.0):
    """Layers of equal optical thickness: d_lo = period n_hi/(n_lo+n_hi), d_hi = period n_lo/(n_lo+n_hi).
    Returns (d_lo, d_hi, f0) with f0 = 1/(4 * optical thickness), the quarter-wave (mid-gap) frequency."""
    d_lo = period * n_hi / (n_lo + n_hi)
    d_hi = period * n_lo / (n_lo + n_hi)
    opt = n_lo * d_lo
    return d_lo, d_hi, 1.0 / (4.0 * opt)


def bloch_cos(f, n_a, d_a, n_b, d_b):
    """cos(K Lambda) for the infinite bilayer: cos da cos db - (1/2)(na/nb + nb/na) sin da sin db."""
    f = np.asarray(f, float)
    da = n_a * TWO_PI * f * d_a
    db = n_b * TWO_PI * f * d_b
    return np.cos(da) * np.cos(db) - 0.5 * (n_a / n_b + n_b / n_a) * np.sin(da) * np.sin(db)


def quarter_wave_gap_edges(n_lo, n_hi, period=1.0):
    """Analytic first stop-band edges of the quarter-wave stack: f0 (1 -+ (2/pi) asin((nh-nl)/(nh+nl)))."""
    _, _, f0 = quarter_wave_bilayer(n_lo, n_hi, period)
    s = (2.0 / np.pi) * np.arcsin((n_hi - n_lo) / (n_hi + n_lo))
    return f0 * (1 - s), f0 * (1 + s)


def quarter_wave_peak_reflectance(n_lo, n_hi, N, n_in=1.0, n_out=1.0):
    """Mid-gap power reflectance of N quarter-wave (lo, hi) pairs, lo adjacent to n_in and hi adjacent to n_out.

    Each quarter-wave layer maps the admittance Y -> n^2 / Y, so after the stack Y = n_out (n_lo/n_hi)^{2N} ... applied
    from the output side gives R = ((n_in n_hi^{2N} - n_out n_lo^{2N}) / (n_in n_hi^{2N} + n_out n_lo^{2N}))^2.
    """
    a = n_in * n_hi ** (2 * N)
    b = n_out * n_lo ** (2 * N)
    return ((a - b) / (a + b)) ** 2


def decay_rate_per_period(f, n_a, d_a, n_b, d_b, period):
    """Power decay in a stop band: T ~ exp(-2 * acosh|cos K Lambda|) per period (nan inside a pass band)."""
    c = np.abs(bloch_cos(f, n_a, d_a, n_b, d_b))
    out = np.full_like(c, np.nan)
    m = c > 1
    out[m] = 2.0 * np.arccosh(c[m])
    return out
