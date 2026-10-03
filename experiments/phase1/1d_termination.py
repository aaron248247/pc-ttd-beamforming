"""Phase 1D: is the narrow in-gap transmission peak a termination (surface) effect? Move the crystal faces and watch it."""
from dataclasses import replace
import numpy as np
from _common import *  # noqa
from pcttd import fdtd2d as F2, io as pio

BASE = F2.Cfg2D(N=8, ny=32, min_time=20000.0, decay_tol=1e-8)
OFFSETS = [-0.15, 0.0, 0.15, 0.30]
import sys
NS = [int(a) for a in sys.argv[1:]] or [8, 10]
OUTNAME = "1d_termination.json" if NS == [8, 10] else "1d_termination_N" + "_".join(map(str, NS)) + ".json"

if __name__ == "__main__":
    cfgs = [replace(BASE, N=n, face_offset=o) for n in NS for o in OFFSETS]
    rs = pool_map(F2.run_2d, cfgs)
    out = dict(config=cfg_dict(BASE), offsets=OFFSETS, Ns=NS, runs=[])
    for c, r in zip(cfgs, rs):
        out["runs"].append(dict(N=c.N, face_offset=c.face_offset, f=r["f"], T=r["T"], R=r["R"], t_end=r["t_end"], hit_max_time=r["hit_max_time"]))
    out["provenance"] = pio.provenance()
    pio.save_json(RES / OUTNAME, out, with_provenance=False)
    print("termination study done")
