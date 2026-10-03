"""Phase 1D: transmission of N-column sections at several resolutions with a long minimum run time (band-edge ring-down).

Saves the full T(f), R(f) curves; edges and attenuation per column are derived from these by 1d_edges.py (no simulation there).
"""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from pcttd import fdtd2d as F2, io as pio

NS = [6, 8, 10, 12, 16]
NYS = [24, 32, 40, 48]
BASE = F2.Cfg2D(N=8, ny=40, min_time=20000.0, decay_tol=1e-8)

if __name__ == "__main__":
    cfgs = [replace(BASE, N=n, ny=ny) for ny in NYS for n in NS]
    rs = pool_map(F2.run_2d, cfgs)
    out = dict(config=cfg_dict(BASE), Ns=NS, nys=NYS, runs=[])
    for c, r in zip(cfgs, rs):
        out["runs"].append(dict(N=c.N, ny=c.ny, res=c.res, f=r["f"], T=r["T"], R=r["R"], t_end=r["t_end"], hit_max_time=r["hit_max_time"], dt=r["dt"]))
    out["provenance"] = pio.provenance()
    pio.save_json(RES / "1d_attenuation.json", out, with_provenance=False)
    print("1D attenuation done")
