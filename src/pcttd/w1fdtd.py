"""Finite W1 line-defect device in MEEP; the plane fields are saved and projected onto guided modes in ``guided``.

Geometry (a = 1, c = 1, TE: Hz out of plane, triangular lattice of air holes in eps = 12; the selected Phase 0.5 family)
    * W1 axis along x; rows j = -J..J (j != 0) at y = j*sqrt(3)/2, row j offset by (j mod 2)/2 in x (see ``geometry.triangular_w1``).
    * Lead region: plain W1 (radius r, rows |j| = 1 unmodified).  Device section: for holes with 0 <= x < N the rows |j| = 1 get
      radius f1*r and are shifted outward by s.  Rows j = +-1 sit at x = p + 1/2, so the section holds exactly N holes per row and
      the structure is mirror symmetric about x = N/2.  The reference planes x = -d_mon and x = N + d_mon are symmetric about it.
    * The mode-analysis cross-section is exactly one MPB supercell (height J*sqrt(3), rows +-J on its edges); the crystal continues
      J_pml more rows inside the y-PML, which only absorbs the evanescent tail.

Why not ``get_eigenmode_coefficients``: in this MEEP build the cross-section mode solve for a periodic (x) structure did not reproduce
the MPB supercell bands (checked: different band frequencies at the same k), whereas the FDTD Bloch phase per period matched MPB.
The decomposition is therefore done in ``guided.project`` with MPB (``w1modes``) mode profiles and the Lorentz-reciprocity overlap.

Pixel alignment: the cell is placed so that pixel centers fall on integer x, hence the reference planes (integer x) coincide with
pixel centers and the periodic mode profile u(0, y) applies to every plane identically.
"""
from __future__ import annotations

import contextlib
import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np

import meep as mp

SQRT3 = math.sqrt(3.0)


@contextlib.contextmanager
def quiet():
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


@dataclass(frozen=True)
class W1Cfg:
    N: int = 8                      # periods in the device section
    r: float = 0.30
    f1: float = 1.0                 # radius factor of rows |j| = 1 inside the device section
    s: float = 0.0                  # outward shift of rows |j| = 1 inside the device section
    eps: float = 12.0
    J: int = 6                      # hole rows on each side of the defect; the analyzed cross-section height is J*sqrt(3)
    J_pml: int = 2                  # extra hole rows continued inside the y-PML
    res: int = 30                   # pixels per lattice constant (integer so integer x lie on pixel centers)
    d_mon: int = 3                  # reference planes this far outside the device section
    d_src: int = 6                  # source plane this far left of the device section
    lead_pad: int = 2               # extra lead between the outermost plane and the x-PML
    dpml_x: int = 3
    dpml_y: float = 1.5
    fcen: float = 0.25
    fwidth: float = 0.40
    freqs: tuple = field(default_factory=lambda: tuple(np.round(np.arange(0.22, 0.275 + 1e-9, 0.001), 6)))
    decay_tol: float = 1e-7
    min_after: float = 150.0
    max_time: float = 30000.0


def hole_list(cfg: W1Cfg):
    """All air holes as (x, y, radius) covering the cell, with the device-section rows modified."""
    x_lo = -(cfg.d_src + cfg.lead_pad + cfg.dpml_x) - 1.0
    x_hi = cfg.N + cfg.d_mon + cfg.lead_pad + cfg.dpml_x + 1.0
    holes = []
    Jt = cfg.J + cfg.J_pml
    for j in list(range(-Jt, 0)) + list(range(1, Jt + 1)):
        off = 0.5 * (j % 2)
        y0 = j * SQRT3 / 2.0
        for p in range(int(math.floor(x_lo)) - 1, int(math.ceil(x_hi)) + 2):
            x = p + off
            if not (x_lo <= x <= x_hi):
                continue
            rad, y = cfg.r, y0
            if abs(j) == 1 and 0.0 <= x < cfg.N:
                rad = cfg.f1 * cfg.r
                y = y0 + math.copysign(1.0, j) * cfg.s        # same convention as geometry.triangular_w1: positive s = outward
            holes.append((x, y, rad))
    return holes


def H_cell(cfg: W1Cfg) -> float:
    return cfg.J * SQRT3


def build(cfg: W1Cfg, src_profile=None):
    """Return the Simulation.  ``src_profile = (y, values)`` shapes the Hz line source across the cross-section (real values)."""
    res = cfg.res
    Lleft = cfg.d_src + cfg.lead_pad + cfg.dpml_x
    Lright = cfg.d_mon + cfg.lead_pad + cfg.dpml_x
    nx = (Lleft + cfg.N + Lright) * res
    xmin = -Lleft - 0.5 / res                        # pixel centers at -Lleft + i/res: integers sit on pixel centers
    sx = nx / res
    ny = int(math.ceil((H_cell(cfg) + 2 * cfg.dpml_y) * res - 1e-9))
    ny += ny % 2                                      # even: y = 0 is a pixel edge, symmetric pixels about the axis
    sy = ny / res
    air = mp.Medium(epsilon=1.0)
    geom = [mp.Cylinder(rad, center=mp.Vector3(x, y), material=air) for x, y, rad in hole_list(cfg)]
    src = mp.GaussianSource(frequency=cfg.fcen, fwidth=cfg.fwidth)
    if src_profile is None:
        amp = lambda p: 1.0 if abs(p.y) < 0.3 else 0.0
    else:
        yy, vv = src_profile
        amp = lambda p: float(np.interp(p.y, yy, vv, left=0.0, right=0.0))
    srcs = [mp.Source(src, component=mp.Hz, center=mp.Vector3(-cfg.d_src, 0), size=mp.Vector3(0, H_cell(cfg)), amp_func=amp)]
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry_center=mp.Vector3(xmin + 0.5 * sx, 0), resolution=res,
                        boundary_layers=[mp.PML(cfg.dpml_x, direction=mp.X), mp.PML(cfg.dpml_y, direction=mp.Y)],
                        geometry=geom, sources=srcs, default_material=mp.Medium(epsilon=cfg.eps), eps_averaging=True)
    return sim


def run_w1(cfg: W1Cfg, src_profile=None) -> dict:
    """Run once; return the DFT plane fields (Ey, Hz) at the input and output reference planes for every requested frequency."""
    sim = build(cfg, src_profile)
    fl = list(cfg.freqs)
    H = H_cell(cfg)
    planes = {}
    for name, x0 in (("in", -cfg.d_mon), ("out", cfg.N + cfg.d_mon)):
        planes[name] = sim.add_dft_fields([mp.Ey, mp.Hz], fl, center=mp.Vector3(x0, 0), size=mp.Vector3(0, H), yee_grid=False)
    meta = {}
    with quiet():
        sim.init_sim()
    for name, x0 in (("in", -cfg.d_mon), ("out", cfg.N + cfg.d_mon)):
        xs, ys, zs, wts = sim.get_array_metadata(dft_cell=planes[name])
        xs = np.atleast_1d(xs)
        if xs.size != 1 or abs(float(xs[0]) - x0) > 1e-9:
            raise RuntimeError(f"monitor plane {name} sampled at x={xs} instead of {x0}")
        meta[name] = (np.asarray(ys, float), np.asarray(wts, float).reshape(-1))
    transit = (cfg.N + cfg.d_src + cfg.d_mon) * 4.0
    with quiet():
        sim.run(until_after_sources=mp.stop_when_dft_decayed(cfg.decay_tol, max(cfg.min_after, transit), cfg.max_time))
    out = dict(f=np.array(fl), t_end=float(sim.meep_time()), dt=float(sim.fields.dt),
               hit_max_time=bool(sim.meep_time() >= cfg.max_time - 1e-9))
    for name in ("in", "out"):
        d = planes[name]
        out[name + "_y"], out[name + "_w"] = meta[name]
        for comp, key in ((mp.Ey, "Ey"), (mp.Hz, "Hz")):
            out[f"{name}_{key}"] = np.array([np.asarray(sim.get_dft_array(d, comp, i)).reshape(-1) for i in range(len(fl))])
    return out
