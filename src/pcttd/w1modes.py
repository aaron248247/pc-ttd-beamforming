"""MPB side of the W1 guided-mode analysis (Phase 2+): even-class dispersion with group index, and guided-mode field profiles.

Units a = 1, c = 1, f = a/lambda. The supercell is the Phase 0.5 one (height M*sqrt(3), rows j = -M..M-1, defect row j = 0 removed).

Parity: MPB's ``run_te_yodd`` is the class with Hz *even* under y -> -y (verified in Phase 2 against ``mpb_bands.defect_bands``
parity = +1: identical frequencies at the same k). This is the class the y-symmetric FDTD source excites.

Folding. The guided band of the plain W1 *decreases* in frequency over k in [0, 0.5] (v_g < 0 there): the physical forward-propagating
wave of that band has Bloch wavenumber k_fwd = -k (equivalently 1 - k). All "forward" quantities below use k_fwd = -kappa, kappa in [0, 0.5],
and the group index is n_g = c/|v_g|.
"""
from __future__ import annotations

import math
import os
import sys
import contextlib

import numpy as np

import meep as mp
from meep import mpb

from . import geometry as G

EPS_HOST = 12.0
AIR = mp.Medium(epsilon=1.0)


@contextlib.contextmanager
def _quiet():
    sys.stdout.flush()
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


def _solver(M, r, f1, s, ks, nb, res, tol, eps=EPS_HOST):
    H, sc = G.triangular_w1(M, r, f1, s)
    geom = [mp.Cylinder(c.radius, material=AIR, center=mp.Vector3(c.x, c.y)) for c in sc]
    ms = mpb.ModeSolver(num_bands=nb, k_points=[mp.Vector3(k) for k in ks], geometry_lattice=mp.Lattice(size=mp.Vector3(1, H)),
                        geometry=geom, default_material=mp.Medium(epsilon=eps), resolution=res, tolerance=tol)
    return ms, H


def even_bands(M=6, r=0.30, f1=1.0, s=0.0, res=32, nk=101, nb=10, tol=1e-9):
    """Hz-even bands on k in [0, 0.5]: freq[nk, nb] and signed v_g[nk, nb] (units of c, Hellmann-Feynman)."""
    ks = np.linspace(0.0, 0.5, nk)
    ms, H = _solver(M, r, f1, s, ks, nb, res, tol)
    freq = np.zeros((nk, nb)); vg = np.zeros((nk, nb))
    st = {"ik": -1}

    def cb(ms_, b):
        if b == 1:
            st["ik"] += 1
            st["vg"] = np.array([v.x for v in ms_.compute_group_velocities()])
            st["f"] = np.array(ms_.get_freqs())
        freq[st["ik"], b - 1] = st["f"][b - 1]
        vg[st["ik"], b - 1] = st["vg"][b - 1]

    with _quiet():
        ms.run_te_yodd(cb)
    return dict(k=ks, freq=freq, vg=vg, H=H, M=M, r=r, f1=f1, s=s, res=res, nk=nk, nb=nb, tol=tol)


def forward_profiles(kappas, band, M=6, r=0.30, f1=1.0, s=0.0, res=32, nb=10, tol=1e-9):
    """Forward guided-mode profiles at the plane x = 0 for k_fwd = -kappa (vg > 0 for a decreasing band).

    Returns dict(y, Ey[nk, ny], Hz[nk, ny], freq[nk]). Fields as given by MPB at the plane x = 0 (Bloch factor 1 there),
    i.e. the periodic part u(0, y); callers add the Bloch factor exp(i k x0) themselves when needed.
    """
    kappas = np.asarray(kappas, float)
    ms, H = _solver(M, r, f1, s, -kappas, max(nb, band), res, tol)
    out = {"Ey": [], "Hz": [], "freq": [], "vg": []}
    st = {}

    def cb(ms_, b):
        if b == 1:
            st["vg"] = np.array([v.x for v in ms_.compute_group_velocities()])
            st["f"] = np.array(ms_.get_freqs())
        if b == band:
            Hf = np.array(ms_.get_hfield(b, bloch_phase=True))[:, :, 0, :]
            Ef = np.array(ms_.get_efield(b, bloch_phase=True))[:, :, 0, :]
            nx = Hf.shape[0]
            i0 = nx // 2                                    # x = (i/nx - 1/2) -> i = nx/2 is x = 0
            out["Hz"].append(Hf[i0, :, 2].copy()); out["Ey"].append(Ef[i0, :, 1].copy())
            out["freq"].append(st["f"][b - 1]); out["vg"].append(st["vg"][b - 1])

    with _quiet():
        ms.run_te_yodd(cb)
    ny = out["Hz"][0].shape[0]
    y = (np.arange(ny) / ny - 0.5) * H
    return dict(y=y, Ey=np.array(out["Ey"]), Hz=np.array(out["Hz"]), freq=np.array(out["freq"]), vg=np.array(out["vg"]), H=H)


def kappa_of_f(bands: dict, band: int, f, vg_floor: float = 2e-3):
    """Invert a monotone-decreasing guided band f(kappa) (kappa in [0, 0.5]) to kappa(f) by a cubic spline through the band samples
    with |v_g| > vg_floor (the flat band ends are excluded).  ``band`` is the 1-based band index within the Hz-even class.
    f must lie inside the retained frequency range."""
    from scipy.interpolate import CubicSpline
    keep = bands["vg"][:, band - 1] < -vg_floor
    fr = bands["freq"][keep, band - 1]; kk = bands["k"][keep]
    if not np.all(np.diff(fr) < 0):
        raise ValueError("band is not monotone decreasing in kappa where |v_g| > vg_floor")
    f = np.asarray(f, float)
    if f.min() < fr.min() - 1e-12 or f.max() > fr.max() + 1e-12:
        raise ValueError("frequency outside the retained band range")
    return CubicSpline(fr[::-1], kk[::-1])(f)
