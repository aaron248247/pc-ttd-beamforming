"""MPB drivers: bulk band structures and W1 supercell dispersion with parity / group velocity / confinement.

Everything dimensionless: a = 1, c = 1, w = a/lambda.

Polarization (MPB and MEEP agree): TE = (Ex, Ey, Hz); TM = (Ez, Hx, Hy).
``Phi`` denotes the out-of-plane field (Hz for TE, Ez for TM); parity is that of Phi under y -> -y about the
defect centre line.
"""
from __future__ import annotations

import contextlib
import math
import os
import sys
import time

import numpy as np

import meep as mp
from meep import mpb

from . import bands as B
from . import geometry as G

EPS_HOST = 12.0
AIR = mp.Medium(epsilon=1.0)

LATTICE_POL = {"triangular": "te", "square": "tm"}


@contextlib.contextmanager
def quiet():
    """Silence C-level stdout (MPB prints per-iteration data regardless of ``verbose``)."""
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


def _run(ms, pol, *fns):
    with quiet():
        (ms.run_te if pol == "te" else ms.run_tm)(*fns)


# --------------------------------------------------------------------------- bulk
def _bulk_solver(lattice: str, r: float, eps: float, res: int, nbands: int, kpts):
    host = mp.Medium(epsilon=eps)
    if lattice == "triangular":
        gl = mp.Lattice(size=mp.Vector3(1, 1), basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
                        basis2=mp.Vector3(math.sqrt(3) / 2, -0.5))
        geom = [mp.Cylinder(r, material=AIR)]
        default = host
    elif lattice == "square":
        gl = mp.Lattice(size=mp.Vector3(1, 1))
        geom = [mp.Cylinder(r, material=host)]
        default = AIR
    else:
        raise ValueError(lattice)
    return mpb.ModeSolver(num_bands=nbands, k_points=kpts, geometry_lattice=gl, geometry=geom,
                          default_material=default, resolution=res)


def high_symmetry_path(lattice: str, per_segment: int = 16):
    V = mp.Vector3
    if lattice == "triangular":     # reciprocal-basis coordinates, MPB tutorial convention
        pts = [V(0, 0), V(0, 0.5), V(-1 / 3, 1 / 3), V(0, 0)]
    else:
        pts = [V(0, 0), V(0.5, 0), V(0.5, 0.5), V(0, 0)]
    return mp.interpolate(per_segment, pts)


def bulk_bands(lattice: str, r: float, eps: float = EPS_HOST, res: int = 32, nbands: int = 6,
               path_per_segment: int = 16, grid_n: int = 14):
    """Bulk bands on the high-symmetry path and on an interior grid covering the whole BZ.

    Returns dict(freqs_path, freqs_grid, ...). The gap is evaluated on the union of both samples.
    """
    pol = LATTICE_POL[lattice]
    path = high_symmetry_path(lattice, path_per_segment)
    grid = [mp.Vector3(i / grid_n, j / grid_n) for i in range(grid_n) for j in range(grid_n)]
    out = {}
    for name, kp in (("path", path), ("grid", grid)):
        ms = _bulk_solver(lattice, r, eps, res, nbands, kp)
        _run(ms, pol)
        out[name] = np.array(ms.all_freqs)
    return out


def bulk_gaps(lattice: str, r: float, eps: float = EPS_HOST, res: int = 32, n_max: int = 3, **kw):
    b = bulk_bands(lattice, r, eps, res, nbands=n_max + 2, **kw)
    both = np.vstack([b["path"], b["grid"]])
    return B.find_gaps(both, n_max=n_max), b


# --------------------------------------------------------------------------- defect supercell
def _geometry(lattice: str, scatterers, eps: float):
    host = mp.Medium(epsilon=eps)
    if lattice == "triangular":
        mat, default = AIR, host
    else:
        mat, default = host, AIR
    geom = [mp.Cylinder(s.radius, material=mat, center=mp.Vector3(s.x, s.y)) for s in scatterers]
    return geom, default


def defect_bands(lattice: str, M: int, r: float, f1: float = 1.0, s: float = 0.0, eps: float = EPS_HOST,
                 res: int = 24, nk: int = 41, gap_n: int = 1, extra: int = 8, tolerance: float = 1e-8,
                 _override=None):
    """Solve a W1 supercell for k in [0, 0.5] and collect per-band data.

    Returns dict with arrays shaped [nk, nbands]: freq, vg (x-component, units of c), parity (+1/-1 of Phi),
    strip and edge energy fractions; plus geometry info and wall time.
    """
    pol = LATTICE_POL[lattice]
    if _override is not None:      # (H, scatterer list): used by verification tests with analytic answers
        H, sc = _override
    else:
        H, sc = (G.triangular_w1 if lattice == "triangular" else G.square_w1)(M, r, f1, s)
    pitch = G.row_pitch(lattice)
    cells = len(sc) + 1
    nb = gap_n * cells + extra
    geom, default = _geometry(lattice, sc, eps)
    ks = np.linspace(0.0, 0.5, nk)
    ms = mpb.ModeSolver(num_bands=nb, k_points=[mp.Vector3(k) for k in ks],
                        geometry_lattice=mp.Lattice(size=mp.Vector3(1, H)), geometry=geom,
                        default_material=default, resolution=res, tolerance=tolerance)

    freq = np.zeros((nk, nb)); vg = np.zeros((nk, nb)); parity = np.zeros((nk, nb), int)
    strip = np.zeros((nk, nb)); edge = np.zeros((nk, nb))
    state = {"ik": -1, "eps": None, "ymask": None}

    def cb(ms_, b):
        if b == 1:
            state["ik"] += 1
            gv = ms_.compute_group_velocities()
            state["vg"] = np.array([v.x for v in gv])
            state["freqs"] = np.array(ms_.get_freqs())
            if state["eps"] is None:
                e = np.array(ms_.get_epsilon())
                nx, ny = e.shape
                y = (np.arange(ny) / ny - 0.5) * H
                state["eps"] = e
                state["strip_mask"] = np.abs(y) <= 1.5 * pitch
                state["edge_mask"] = np.abs(y) >= H / 2 - pitch
                state["flip"] = (ny - np.arange(ny)) % ny
                # the structure itself must be mirror symmetric on the grid (checked loosely: subpixel smoothing)
                state["asym"] = float(np.max(np.abs(e - e[:, state["flip"]])))
        ik = state["ik"]
        freq[ik, b - 1] = state["freqs"][b - 1]
        vg[ik, b - 1] = state["vg"][b - 1]
        Hf = np.array(ms_.get_hfield(b))[:, :, 0, :]
        Df = np.array(ms_.get_dfield(b))[:, :, 0, :]
        phi = Hf[:, :, 2] if pol == "te" else np.array(ms_.get_efield(b))[:, :, 0, 2]
        nrm = np.vdot(phi, phi).real
        if nrm < 1e-20:                      # e.g. the w=0 band at k=0: parity undefined, fixed up below
            parity[ik, b - 1] = 0
        else:
            ov = np.vdot(phi, phi[:, state["flip"]]) / nrm
            parity[ik, b - 1] = +1 if ov.real >= 0 else -1
        u = (np.abs(Df) ** 2).sum(-1) / state["eps"] + (np.abs(Hf) ** 2).sum(-1)
        tot = u.sum()
        strip[ik, b - 1] = u[:, state["strip_mask"]].sum() / tot
        edge[ik, b - 1] = u[:, state["edge_mask"]].sum() / tot

    t0 = time.time()
    _run(ms, pol, cb)
    wall = time.time() - t0
    # undefined-parity entries (vanishing Phi) inherit the parity of the same band at the nearest k with a field
    for ik, b in zip(*np.nonzero(parity == 0)):
        good = np.nonzero(parity[:, b] != 0)[0]
        parity[ik, b] = parity[good[np.argmin(np.abs(good - ik))], b]
    return dict(k=ks, freq=freq, vg=vg, parity=parity, strip=strip, edge=edge, H=H, n_scatterers=len(sc),
                r_min=min((s_.radius for s_ in sc), default=float('nan')), wall_s=wall, eps_mirror_asym=state["asym"],
                lattice=lattice, M=M, r=r, f1=f1, s=s, eps=eps, res=res, nk=nk)


# --------------------------------------------------------------------------- per-design analysis
def analyze_design(d: dict, gap: B.Gap, n_levels: int = 1500) -> dict:
    """Apply the pre-declared Phase 0.5 definitions to one defect-band calculation."""
    k = d["k"]
    br = B.class_branches(d["freq"], d["vg"], d["parity"], d["strip"], d["edge"])
    levels = np.linspace(gap.w_lo, gap.w_hi, n_levels + 2)[1:-1]
    c = B.find_crossings(k, br, levels)
    sa, se, eo = B.single_mode_masks(c)
    cnt = c.counts()

    res = dict(gap_w_lo=gap.w_lo, gap_w_hi=gap.w_hi, gap_rel=gap.rel, max_modes_in_gap=int(cnt.max()),
               n_defect_branches_in_gap=int(sum(((b["w"] > gap.w_lo) & (b["w"] < gap.w_hi)).any() for b in br)))

    def longest(mask):
        rr = B.runs(mask)
        if not rr:
            return None, 0.0
        fw = [B.frac_width(levels, a, b) for a, b in rr]
        i = int(np.argmax(fw))
        return rr[i], fw[i]

    run_se, fw_se = longest(se)
    _, fw_sa = longest(sa)
    _, fw_eo = longest(eo)
    res.update(sm_strict_even_frac=fw_se, sm_strict_any_frac=fw_sa, sm_even_only_frac=fw_eo)

    # velocity/strip/edge of the single even forward crossing, on levels where strict_even
    vg_l, strip_l, edge_l = B.per_level_single(c, se[c.level_idx])
    ng_l = np.where(np.isfinite(vg_l) & (vg_l != 0), 1.0 / np.abs(vg_l), np.nan)   # n_g = c/|v_g|
    res["ng_max_strict_even"] = float(np.nanmax(ng_l)) if np.isfinite(ng_l).any() else 0.0

    gbps = {}
    for n_min in (5, 10, 20):
        best = None
        for a, b_ in B.runs(se):
            if b_ - a < 3:
                continue
            r_ = B.best_flat_window(levels[a:b_], ng_l[a:b_], n_min)
            if r_ is not None and (best is None or r_["gbp"] > best["gbp"]):
                best = r_
        gbps[str(n_min)] = best
    res["flat_window"] = gbps
    res["gbp10"] = gbps["10"]["gbp"] if gbps["10"] else 0.0

    if run_se is not None:
        a, b_ = run_se
        mid = (a + b_ - 1) // 2
        res["window"] = dict(w_lo=float(levels[a]), w_hi=float(levels[b_ - 1]), w_center=float(levels[mid]),
                             ng_at_center=float(ng_l[mid]), ng_min=float(np.nanmin(ng_l[a:b_])),
                             ng_max=float(np.nanmax(ng_l[a:b_])), strip_at_center=float(strip_l[mid]),
                             edge_max=float(np.nanmax(edge_l[a:b_])),
                             vg_sign_constant=bool(np.all(np.sign(vg_l[a:b_]) == np.sign(vg_l[a]))),
                             vg_sign=int(np.sign(vg_l[a])))
    else:
        res["window"] = None

    # C6: even strict window >= 1% wide that stays >= 5% of the gap width from both gap edges
    margin = 0.05 * gap.width
    m = se & (levels >= gap.w_lo + margin) & (levels <= gap.w_hi - margin)
    _, fw_m = longest(m)
    res["sm_strict_even_frac_margin"] = fw_m
    res["c6_pass"] = bool(fw_m >= 0.01)

    # diagnostic: MPB (Hellmann-Feynman) v_g vs finite-difference slope of w(k), in-gap points only
    errs = []
    for b in br:
        fd = np.gradient(b["w"], k)
        ing = (b["w"] > gap.w_lo) & (b["w"] < gap.w_hi) & (np.abs(b["vg"]) > 0.02)
        ing[[0, -1]] = False
        errs.extend(list(np.abs(fd[ing] - b["vg"][ing]) / np.abs(b["vg"][ing])))
    res["vg_fd_rel_err"] = dict(n=len(errs), median=float(np.median(errs)) if errs else None,
                                p95=float(np.percentile(errs, 95)) if errs else None,
                                max=float(np.max(errs)) if errs else None)
    return res
