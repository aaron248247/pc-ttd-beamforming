"""One-dimensional MEEP driver (fields Ex, Hy along z) for propagation and multilayer transmission experiments.

All lengths in units of a (c = 1). Positions are snapped to the Yee grid so monitor separations are exact.
The DFT frequencies are passed as an explicit list and read back from the DFT object, so the frequency axis used by
the extraction code is the one MEEP actually used (an earlier probe that assumed MEEP's (fcen, df, nfreq) arguments
were (fmin, fmax, nfreq) produced delays off by exactly the grid-spacing ratio; see tests/test_fdtd_1d.py).
"""
from __future__ import annotations

import contextlib
import os
import sys
from dataclasses import dataclass, field, replace

import numpy as np

import meep as mp


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
class Cfg1D:
    res: float = 40.0
    dpml: float = 2.0
    z_lo: float = -8.0                    # interior (non-PML) extent
    z_hi: float = 18.0
    z_src: float = -5.0
    monitors: tuple = (0.0, 10.0)
    freqs: tuple = tuple(np.round(np.arange(0.15, 0.45 + 1e-9, 0.002), 6))
    src_fcen: float = 0.30
    src_fwidth: float = 0.30
    blocks: tuple = ()                    # ((z_left, z_right, eps), ...)
    eps_bg: float = 1.0
    decay_tol: float = 1e-10
    min_time: float | None = None         # None -> automatic: transit time to the farthest monitor + 20 (see run_1d)
    max_time: float = 20000.0
    decimation: int = 1
    yee_grid: bool = False       # True returns zeros for point monitors; interpolated value has the same phase (real cos(k dx/2) factor)

    def snap(self, z):
        return round(z * self.res) / self.res


def run_1d(cfg: Cfg1D) -> dict:
    """Run one 1-D simulation; return complex DFT of Ex at each monitor for every frequency."""
    dx = 1.0 / cfg.res
    length = (cfg.z_hi - cfg.z_lo) + 2 * cfg.dpml
    zc = 0.5 * (cfg.z_lo + cfg.z_hi)
    geom = [mp.Block(size=mp.Vector3(mp.inf, mp.inf, zr - zl), center=mp.Vector3(0, 0, 0.5 * (zl + zr)),
                     material=mp.Medium(epsilon=eps)) for zl, zr, eps in cfg.blocks]
    src = mp.GaussianSource(frequency=cfg.src_fcen, fwidth=cfg.src_fwidth)
    sim = mp.Simulation(cell_size=mp.Vector3(0, 0, length), geometry_center=mp.Vector3(0, 0, zc), dimensions=1,
                        resolution=cfg.res, boundary_layers=[mp.PML(cfg.dpml)] if cfg.dpml > 0 else [],
                        sources=[mp.Source(src, component=mp.Ex, center=mp.Vector3(0, 0, cfg.snap(cfg.z_src)))],
                        geometry=geom, default_material=mp.Medium(epsilon=cfg.eps_bg))
    freqs = list(cfg.freqs)
    objs = [sim.add_dft_fields([mp.Ex], freqs, decimation_factor=cfg.decimation, yee_grid=cfg.yee_grid,
                               where=mp.Volume(center=mp.Vector3(0, 0, cfg.snap(z)), size=mp.Vector3()))
            for z in cfg.monitors]
    # A decay test on a DFT that has not started yet passes trivially, so never stop before the wave can have arrived:
    # transit time at the slowest speed present (c / n_max) to the farthest monitor, plus a margin for the pulse.
    n_max = float(np.sqrt(max([cfg.eps_bg] + [b[2] for b in cfg.blocks])))
    t_transit = n_max * max(abs(z - cfg.z_src) for z in cfg.monitors)
    min_time = cfg.min_time if cfg.min_time is not None else t_transit + 20.0
    with _quiet():
        sim.run(until_after_sources=mp.stop_when_dft_decayed(cfg.decay_tol, min_time, cfg.max_time))
    E = np.array([[complex(sim.get_dft_array(o, mp.Ex, i)) for i in range(len(freqs))] for o in objs])
    f_used = np.array(list(objs[0].freq))
    return dict(f=f_used, E=E, dt=float(sim.fields.dt), t_end=float(sim.meep_time()), dx=dx,
                hit_max_time=bool(sim.meep_time() >= cfg.max_time - 1e-9), min_time=float(min_time), t_transit=float(t_transit))


def delay_metrics(f, H, tau_ideal, tau_pred=None, band=None):
    """Compare extracted group delays to an ideal (and optionally a Yee-predicted) value.

    ``tau_ideal`` / ``tau_pred`` may be scalars or arrays over ``f``. Both extraction implementations are reported.
    """
    from . import extract as X
    f = np.asarray(f, float)
    tp, ok = X.group_delay_phase(f, H)
    tr = X.group_delay_ratio(f, H)
    m = np.ones_like(f, bool) if band is None else X.band_mask(f, *band)
    m &= np.isfinite(tp)
    ideal = np.broadcast_to(np.asarray(tau_ideal, float), f.shape)
    out = dict(unwrap_ok=bool(ok), n_points=int(m.sum()),
               mean_err=float(np.mean(tp[m] - ideal[m])), max_abs_err=float(np.max(np.abs(tp[m] - ideal[m]))),
               max_abs_err_ratio_method=float(np.max(np.abs(tr[m][1:-1] - ideal[m][1:-1]))),
               max_phase_vs_ratio=float(np.max(np.abs(tp[m][1:-1] - tr[m][1:-1]))),
               mean_tau=float(np.mean(tp[m])))
    if tau_pred is not None:
        pred = np.broadcast_to(np.asarray(tau_pred, float), f.shape)
        out["max_abs_resid_vs_yee"] = float(np.max(np.abs(tp[m] - pred[m])))
        out["mean_resid_vs_yee"] = float(np.mean(tp[m] - pred[m]))
    return out
