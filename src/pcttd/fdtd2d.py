"""Two-dimensional MEEP driver: normal-incidence plane wave on a finite section of the triangular air-hole crystal.

Geometry (a = 1; all outside material is vacuum). Propagation along +x (the nearest-neighbor, Gamma-K direction).
The cell has height sqrt(3) and is periodic in y, so a y-uniform Hz line source is an infinite plane wave at normal incidence and the
crystal is exactly the rectangular two-row cell that MPB solves at k_y = 0.
Crystal: a block of dielectric host (eps = 12) spanning x in [-0.5, N-0.5] with air holes (radius r) at
    y = 0:           x = 0, 1, ..., N-1
    y = +-sqrt(3)/2: x = 0.5, ..., N-1.5      (both periodic images placed explicitly)
i.e. symmetric about its center and terminated by full y = 0 hole rows; N-1 lattice periods between end-hole centers.
Polarization: TE (Hz out of plane; Ex, Ey in plane), as selected in Phase 0.5.
T and R from flux monitors (explicit frequency list), normalized by an empty (vacuum) run with identical monitors.
"""
from __future__ import annotations

import contextlib
import math
import os
import sys
from dataclasses import dataclass

import numpy as np

import meep as mp

SQRT3 = math.sqrt(3.0)


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


@dataclass(frozen=True)
class Cfg2D:
    N: int = 12
    r: float = 0.30
    eps: float = 12.0
    ny: int = 32                      # pixels across the periodic height sqrt(3): resolution = ny/sqrt(3) (cell height exact)
    dpml: float = 2.0
    src_dist: float = 3.0             # source this far left of the crystal's left face (x = -0.5)
    mon_dist: float = 1.0             # flux monitors this far outside the faces
    pad: float = 1.0                  # vacuum between the PML inner edge and the source / transmission monitor
    freqs: tuple = tuple(np.round(np.arange(0.12, 0.36 + 1e-9, 0.002), 6))
    src_fcen: float = 0.24
    src_fwidth: float = 0.40
    decay_tol: float = 1e-8
    max_time: float = 60000.0
    min_time: float | None = None       # None -> transit time + 30; set larger to outlast band-edge ring-down
    decimation: int = 1
    face_offset: float = 0.0           # extra host material beyond each end face (termination study); faces at -0.5-offset, N-0.5+offset

    @property
    def res(self):
        return self.ny / SQRT3

    def faces(self):
        return -0.5 - self.face_offset, self.N - 0.5 + self.face_offset


def _build(cfg: Cfg2D, with_crystal: bool):
    xl_face, xr_face = cfg.faces()
    x_src = xl_face - cfg.src_dist
    xl = x_src - cfg.pad
    xr = xr_face + cfg.mon_dist + cfg.pad
    res = cfg.ny / SQRT3
    sx = math.ceil(((xr - xl) + 2 * cfg.dpml) * res - 1e-9) / res        # whole number of pixels along x
    cx = 0.5 * (xl + xr)
    geom = []
    if with_crystal:
        geom.append(mp.Block(size=mp.Vector3(xr_face - xl_face, mp.inf), center=mp.Vector3(0.5 * (xl_face + xr_face), 0),
                             material=mp.Medium(epsilon=cfg.eps)))
        air = mp.Medium(epsilon=1.0)
        for p in range(cfg.N):
            geom.append(mp.Cylinder(cfg.r, center=mp.Vector3(float(p), 0.0), material=air))
        for p in range(cfg.N - 1):
            for sy in (+1, -1):
                geom.append(mp.Cylinder(cfg.r, center=mp.Vector3(p + 0.5, sy * SQRT3 / 2), material=air))
    src = mp.GaussianSource(frequency=cfg.src_fcen, fwidth=cfg.src_fwidth)
    sim = mp.Simulation(cell_size=mp.Vector3(sx, SQRT3), geometry_center=mp.Vector3(cx, 0), resolution=cfg.ny / SQRT3,
                        boundary_layers=[mp.PML(cfg.dpml, direction=mp.X)], geometry=geom,
                        sources=[mp.Source(src, component=mp.Hz, center=mp.Vector3(x_src, 0), size=mp.Vector3(0, SQRT3))],
                        k_point=mp.Vector3(0, 0, 0), default_material=mp.Medium(epsilon=1.0))
    return sim


def _flux(sim, cfg, x):
    return sim.add_flux(list(cfg.freqs), mp.FluxRegion(center=mp.Vector3(x, 0), size=mp.Vector3(0, SQRT3)),
                        decimation_factor=cfg.decimation)


def _min_time(cfg):
    # transit through the crystal at the slowest speed present (c / sqrt(eps)) plus margin for the pulse
    xl_face, xr_face = cfg.faces()
    return math.sqrt(cfg.eps) * (cfg.src_dist + (xr_face - xl_face) + cfg.mon_dist) + 30.0


def run_2d(cfg: Cfg2D) -> dict:
    xl_face, xr_face = cfg.faces()
    xr_mon, xt_mon = xl_face - cfg.mon_dist, xr_face + cfg.mon_dist
    mt = cfg.min_time if cfg.min_time is not None else _min_time(cfg)
    stop = lambda: mp.stop_when_dft_decayed(cfg.decay_tol, mt, cfg.max_time)
    # 1. empty run: incident flux and the incident field at the reflection plane
    s0 = _build(cfg, False)
    refl0, tran0 = _flux(s0, cfg, xr_mon), _flux(s0, cfg, xt_mon)
    with _quiet():
        s0.run(until_after_sources=stop())
    inc = np.array(mp.get_fluxes(tran0))
    inc_refl = np.array(mp.get_fluxes(refl0))
    f_used = np.array(mp.get_flux_freqs(tran0))
    data0 = s0.get_flux_data(refl0)
    t0 = s0.meep_time()
    # 2. crystal run
    s1 = _build(cfg, True)
    refl1, tran1 = _flux(s1, cfg, xr_mon), _flux(s1, cfg, xt_mon)
    s1.load_minus_flux_data(refl1, data0)
    with _quiet():
        s1.run(until_after_sources=stop())
    T = np.array(mp.get_fluxes(tran1)) / inc
    R = -np.array(mp.get_fluxes(refl1)) / inc
    return dict(f=f_used, T=T, R=R, inc=inc, inc_refl=inc_refl, t_end_empty=float(t0), t_end=float(s1.meep_time()),
                dt=float(s1.fields.dt),
                hit_max_time=bool(max(t0, s1.meep_time()) >= cfg.max_time - 1e-9))
